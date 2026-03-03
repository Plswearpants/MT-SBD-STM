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
    num_slices    = size(Y, 3);

    % ---------------------------------------------------------------------
    % Determine trusted-slice step weights
    % ---------------------------------------------------------------------
    Y_used = Y;
    A1_used = A1_all_matrix;

    trusted_ratio_threshold = input(sprintf('Enter trusted-slice std-ratio threshold (e.g. %.2f): ', ...
        cfg.blockRun.trusted_ratio_threshold_default));
    if isempty(trusted_ratio_threshold)
        trusted_ratio_threshold = cfg.blockRun.trusted_ratio_threshold_default;
    end

    manual_trusted_slices = cell(1, num_kernels);
    if cfg.blockRun.use_default_manual_trusted_slices && num_kernels >= 5
        manual_trusted_slices{1} = [1,4,5,8,10];
        manual_trusted_slices{2} = [1,5,8,9,10];
        manual_trusted_slices{3} = 7:11;
        manual_trusted_slices{4} = 7:11;
        manual_trusted_slices{5} = [3,5,6,7];
    end

    [params.slice_weights, params.slice_weight_details] = ...
        build_auto_trusted_slice_weights(A1_used, eta_data3d, trusted_ratio_threshold, ...
        cfg.blockRun.show_trusted_plot, manual_trusted_slices);

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
    params.lambda1_unweighted = sqrt(num_slices) .* params.lambda1_base;
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
    params.noise_var       = eta_data3d;
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
    fprintf('Trusted-slice weights ready. Counts per kernel: %s\n', mat2str(trusted_counts));
    fprintf('Lambda weighted (sqrt(trusted_count)): %s\n', mat2str(params.lambda1_weighted, 4));
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
    eta3dall = permute(repmat(eta_data3d, [outerloop_maxIT, 1]), [2,1]);
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

    meta.stage = "run";

end

