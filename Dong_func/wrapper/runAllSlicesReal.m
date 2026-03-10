function [log, data, params, meta, cfg] = runAllSlicesReal(log, data, params, meta, cfg)
%RUNALLSLICESREAL Run MT-SBD on all slices (block run).
%
%   [log, data, params, meta, cfg] = runAllSlicesReal(log, data, params, meta, cfg)
%
%   This wrapper encapsulates the core "block run" logic from
%   scripts/MTSBD_block_realdata1.m:
%       - builds A1_used = A1_all_matrix and Y_used = Y
%       - computes trusted-slice weights via build_auto_trusted_slice_weights
%       - sets lambda1_base and weighted/unweighted variants
%       - configures params for MTSBD_all_slice or MTSBD_Xregulated_all_slices
%       - runs the chosen algorithm
%       - computes observation_fidelity
%       - stores results under data.real.blockRun
%
%   Visualization-heavy and experimental sections from the legacy script
%   (movies, additional padded runs, sequential runs) are intentionally
%   omitted.

    arguments
        log  struct
        data struct
        params struct
        meta struct
        cfg  struct
    end

    if ~isfield(data, "real") || ~isfield(data.real, "proliferation")
        error('runAllSlicesReal: proliferation results missing. Run initProliferationReal first.');
    end
    if ~isfield(data.real, "Y")
        error('runAllSlicesReal: data.real.Y is missing. Run preprocessRealData first.');
    end

    Y = data.real.Y;

    A1_all_matrix = data.real.proliferation.A1_all_matrix;
    eta_data3d    = data.real.proliferation.eta_data3d;
    kernel_sizes  = data.real.ref.kernel_sizes;
    num_kernels   = params.refSlice.num_kernels;
    num_slices_full = size(Y, 3);

    % ---------------------------------------------------------------------
    % Optional: run only a chosen subset of slices
    % ---------------------------------------------------------------------
    slice_indices = 1:num_slices_full;
    if isfield(params, "blockRun") && isfield(params.blockRun, "slices_to_run") ...
            && ~isempty(params.blockRun.slices_to_run)
        slice_indices = params.blockRun.slices_to_run;
    elseif isfield(params, "blockRun") && isfield(params.blockRun, "interactive") ...
            && params.blockRun.interactive
        slice_indices_in = input('Enter slice indices to run (e.g. 1:10 or [1 3 5]; empty = all): ');
        if ~isempty(slice_indices_in)
            slice_indices = slice_indices_in;
        end
    end

    if islogical(slice_indices)
        if numel(slice_indices) ~= num_slices_full
            error('params.blockRun.slices_to_run logical mask must have length equal to number of slices.');
        end
        slice_indices = find(slice_indices);
    end

    if ~isnumeric(slice_indices) || isempty(slice_indices)
        error('params.blockRun.slices_to_run must be a numeric vector of slice indices or a logical mask.');
    end
    slice_indices = slice_indices(:).';
    if any(mod(slice_indices, 1) ~= 0)
        error('params.blockRun.slices_to_run must contain integer indices.');
    end
    if any(slice_indices < 1) || any(slice_indices > num_slices_full)
        error('params.blockRun.slices_to_run indices must be within 1..%d.', num_slices_full);
    end
    slice_indices = unique(slice_indices, 'stable');

    num_slices_run = numel(slice_indices);
    if num_slices_run ~= num_slices_full
        fprintf('Running block algorithm on slice subset (%d/%d): %s\n', ...
            num_slices_run, num_slices_full, mat2str(slice_indices));
    else
        fprintf('Running block algorithm on all slices (%d).\n', num_slices_full);
    end

    % ---------------------------------------------------------------------
    % Determine trusted-slice step weights
    % ---------------------------------------------------------------------
    Y_used = Y(:,:,slice_indices);
    A1_used = A1_all_matrix;
    for k = 1:num_kernels
        A1_used{k} = A1_all_matrix{k}(:,:,slice_indices);
    end
    eta_data3d_used = eta_data3d(slice_indices);

    % Optional: trusted-slice weighting (can be disabled if helper missing)
    use_trusted_weights = true;
    if isfield(cfg, "blockRun") && isfield(cfg.blockRun, "use_trusted_slice_weights") ...
            && ~isempty(cfg.blockRun.use_trusted_slice_weights)
        use_trusted_weights = logical(cfg.blockRun.use_trusted_slice_weights);
    end
    if isfield(params, "blockRun") && isfield(params.blockRun, "use_trusted_slice_weights") ...
            && ~isempty(params.blockRun.use_trusted_slice_weights)
        use_trusted_weights = logical(params.blockRun.use_trusted_slice_weights);
    end

    trusted_ratio_threshold = cfg.blockRun.trusted_ratio_threshold_default;
    if use_trusted_weights
        if isfield(params, "blockRun") && isfield(params.blockRun, "interactive") ...
                && params.blockRun.interactive
            trusted_ratio_threshold_in = input(sprintf('Enter trusted-slice std-ratio threshold (e.g. %.2f): ', ...
                cfg.blockRun.trusted_ratio_threshold_default));
            if ~isempty(trusted_ratio_threshold_in)
                trusted_ratio_threshold = trusted_ratio_threshold_in;
            end
        end

        manual_trusted_slices = cell(1, num_kernels);
        if cfg.blockRun.use_default_manual_trusted_slices && num_kernels >= 5
            manual_trusted_slices{1} = [1,4,5,8,10];
            manual_trusted_slices{2} = [1,5,8,9,10];
            manual_trusted_slices{3} = 7:11;
            manual_trusted_slices{4} = 7:11;
            manual_trusted_slices{5} = [3,5,6,7];
        end

        % manual_trusted_slices values are specified in FULL slice indices; map to
        % positions within the chosen slice subset for build_auto_trusted_slice_weights.
        manual_trusted_slices_used = manual_trusted_slices;
        if num_slices_run ~= num_slices_full
            for k = 1:num_kernels
                if isempty(manual_trusted_slices{k})
                    manual_trusted_slices_used{k} = [];
                    continue;
                end
                abs_idx = manual_trusted_slices{k}(:).';
                [tf, loc] = ismember(abs_idx, slice_indices);
                manual_trusted_slices_used{k} = loc(tf);
            end
        end

        if exist('build_auto_trusted_slice_weights', 'file') ~= 2
            warning(['build_auto_trusted_slice_weights.m not found on path. ', ...
                'Falling back to unweighted mode (all slices treated as trusted).']);
            use_trusted_weights = false;
        else
            [params.slice_weights, params.slice_weight_details] = ...
                build_auto_trusted_slice_weights(A1_used, eta_data3d_used, trusted_ratio_threshold, ...
                cfg.blockRun.show_trusted_plot, manual_trusted_slices_used);
        end
    end

    if ~use_trusted_weights
        % Unweighted: allow MTSBD_all_slice to default slice_weights = ones
        params.slice_weights = [];
        params.slice_weight_details = struct();
        params.slice_weight_details.trusted_counts = num_slices_run * ones(1, num_kernels);
        params.slice_weight_details.method = "unweighted_all_trusted";
        params.slice_weight_details.trusted_ratio_threshold = trusted_ratio_threshold;
    end

    % ---------------------------------------------------------------------
    % Configure block-run parameters
    % ---------------------------------------------------------------------
    miniloop_iteration = cfg.blockRun.miniloop_iteration;
    outerloop_maxIT    = cfg.blockRun.outerloop_maxIT;

    params.lambda1_base = cfg.blockRun.lambda1_base;
    if numel(params.lambda1_base) ~= num_kernels
        error('params.lambda1_base must have one value per kernel.');
    end

    trusted_counts = params.slice_weight_details.trusted_counts;
    params.lambda1_weighted   = sqrt(trusted_counts) .* params.lambda1_base;
    params.lambda1_unweighted = sqrt(num_slices_run) .* params.lambda1_base;
    params.lambda1            = params.lambda1_unweighted;

    params.phase2   = cfg.blockRun.phase2;
    params.kplus    = ceil(cfg.blockRun.kplus_factor * kernel_sizes);
    params.lambda2  = cfg.blockRun.lambda2;
    params.nrefine  = cfg.blockRun.nrefine;
    params.signflip = cfg.blockRun.signflip;
    params.xpos     = cfg.blockRun.xpos;
    params.getbias  = cfg.blockRun.getbias;
    params.Xsolve   = cfg.blockRun.Xsolve;
    params.use_Xregulated = cfg.blockRun.use_Xregulated;
    params.noise_var       = eta_data3d_used;
    params.kernel_update_order = 1:num_kernels;

    use_custom_order = false;
    if cfg.blockRun.allow_custom_update_order
        use_custom_order = input('Use custom kernel update order for MTSBD_all_slice? (0/1): ');
    end
    if ~isempty(use_custom_order) && use_custom_order
        custom_order = input(sprintf('Enter kernel update permutation of 1:%d (e.g. [2 1 3 ...]): ', num_kernels));
        custom_order = custom_order(:).';
        if numel(custom_order) ~= num_kernels || any(custom_order < 1) || ...
                any(custom_order > num_kernels) || numel(unique(custom_order)) ~= num_kernels
            error('Invalid kernel update order. Must be a permutation of 1:num_kernels.');
        end
        params.kernel_update_order = custom_order;
    end
    fprintf('Kernel update order: %s\n', mat2str(params.kernel_update_order));

    kernel_sizes_single = squeeze(max(kernel_sizes, [], 1));
    if use_trusted_weights
        fprintf('Trusted-slice weights ready. Counts per kernel: %s\n', mat2str(trusted_counts));
        fprintf('Lambda weighted (sqrt(trusted_count)): %s\n', mat2str(params.lambda1_weighted, 4));
    else
        fprintf('Trusted-slice weighting disabled (unweighted mode; all slices treated as trusted).\n');
    end
    fprintf('Lambda unweighted (sqrt(total_slices)): %s\n', mat2str(params.lambda1_unweighted, 4));

    % ---------------------------------------------------------------------
    % Set up display functions
    % ---------------------------------------------------------------------
    figure;
    dispfun = cell(1, num_kernels);
    for n = 1:num_kernels
        dispfun{n} = @(Y_, A, X, kernel_sizes_sing, kplus) ... %#ok<NASGU,INUSD>
            showims(Y_used, A1_used{n}, data.real.ref.X_ref(:,:,n), A, X, kernel_sizes_single, kplus, 1);
    end

    % ---------------------------------------------------------------------
    % Run block algorithm
    % ---------------------------------------------------------------------
    if params.use_Xregulated
        [REG_Aout_ALL, REG_Xout_ALL, REG_bout_ALL, REG_extras_ALL] = ...
            MTSBD_Xregulated_all_slices(Y_used, kernel_sizes, params, dispfun, ...
            A1_used, miniloop_iteration, outerloop_maxIT); %#ok<NASGU,INUSD>
        error('X-regulated variant not yet wired into data.real storage. Use non-regulated path for now.');
    else
        [Aout_ALL, Xout_ALL, bout_ALL, ALL_extras] = ...
            MTSBD_all_slice(Y_used, kernel_sizes, params, dispfun, ...
            A1_used, miniloop_iteration, outerloop_maxIT);
    end

    % Observation fidelity
    eta3dall = permute(repmat(eta_data3d_used, [outerloop_maxIT, 1]), [2,1]);
    observation_fidelity = eta3dall ./ squeeze(var(ALL_extras.phase1.residuals, 0, [1,2]));

    % ---------------------------------------------------------------------
    % Store results
    % ---------------------------------------------------------------------
    if ~isfield(data.real, "blockRun")
        data.real.blockRun = struct();
    end

    data.real.blockRun.Aout_ALL           = Aout_ALL;
    data.real.blockRun.Xout_ALL           = Xout_ALL;
    data.real.blockRun.bout_ALL           = bout_ALL;
    data.real.blockRun.ALL_extras         = ALL_extras;
    data.real.blockRun.observation_fidelity = observation_fidelity;
    data.real.blockRun.trusted_counts     = trusted_counts;
    data.real.blockRun.slice_indices      = slice_indices;
    data.real.blockRun.num_slices_full    = num_slices_full;

    meta.stage = "run";

    % ---------------------------------------------------------------------
    % Save block-run checkpoint (optional)
    % ---------------------------------------------------------------------
    if isfield(params, "blockRun") && isfield(params.blockRun, "save_checkpoint") ...
            && params.blockRun.save_checkpoint
        if isfield(cfg, "io") && isfield(cfg.io, "blockrun_output_file") ...
                && ~isempty(cfg.io.blockrun_output_file)
            blockrun_file = cfg.io.blockrun_output_file;
        else
            base = 'realdata';
            if isfield(cfg, "load") && isfield(cfg.load, "data_file") && ~isempty(cfg.load.data_file)
                [~, base, ~] = fileparts(cfg.load.data_file);
            end
            blockrun_file = sprintf('%s_blockrun_checkpoint.mat', base);
        end
        % Avoid overwriting an existing checkpoint: if the target file
        % already exists, append a timestamp suffix.
        if exist(blockrun_file, 'file')
            [fpath, fname, fext] = fileparts(blockrun_file);
            ts = datestr(now, 'yyyymmdd_HHMMSS');
            blockrun_file = fullfile(fpath, sprintf('%s_%s%s', fname, ts, fext));
        end
        save(blockrun_file, 'log', 'data', 'params', 'meta', 'cfg', '-v7.3');
        fprintf('Saved block-run checkpoint to %s.\n', blockrun_file);
    end

end

