function [data, params] = initializeBlock(log, data, params, varargin)
%INITIALIZEBLOCK Unified block initialization: proliferation or block_manual (IK01A-style for block).
%
%   Single entry point for initializing kernels across all slices. Dispatches
%   on params.block.block_init_method: 'proliferation' (IS01A centers) or
%   'block_manual' (ref-slice manual/random centers propagated to all slices).
%   Writes A1_all, A1_all_matrix, init_kernel_centers to data.block (via organizeData).
%
%   [data, params] = initializeBlock(log, data, params, ...)
%
%   INPUTS:
%       log                 - Log struct with .path and .file
%       data                - Data from previous blocks (hierarchical or flat)
%       params              - Params from previous blocks (hierarchical or flat)
%
%       OPTIONAL (Name-Value), override params.block:
%       'block_init_method' - 'proliferation' | 'block_manual' (default: 'proliferation')
%       'A1_matrix_unify_size' - 'max_per_kernel' (default) | 'ref_slice'
%       'window_type_proliferation' - Window for proliferated kernels (default: {'gaussian', 2.5})
%       'interactive_size_adjust'   - Allow interactive kernel size adjustment (default: false)
%       'use_matrix_format'        - Output A1_all_matrix (default: true)
%       'show_reference_kernels'   - Display ref-slice kernels (default: true)
%
%   OUTPUTS:
%       data  - Updated with data.block.A_all_init, A_all_init_matrix, init_kernel_centers
%       params - Updated with params.block.block_init_method and related
%
%   REQUIRED (method = 'proliferation'): data.most_isolated_points, params.target_kernel_sizes (IS01A).
%   REQUIRED (method = 'block_manual'):  Y, Y_ref, params (num_kernels, num_slices, ref_slice, kernel_sizes).
%       Previous data.slice is not used; manual selection is always run and overwrites data.slice.
%
%   See also: initializeProliferation, initializeKernelsRef, initialize_kernels_proliferation, organizeData

    % Validate inputs
    if ~isstruct(log) || ~isfield(log, 'path') || ~isfield(log, 'file')
        error('log must be a struct with .path and .file fields');
    end
    if ~isstruct(data) || ~isstruct(params)
        error('data and params must be structs');
    end

    % Parse optional name-value (override params.block)
    p = inputParser;
    addParameter(p, 'block_init_method', 'proliferation', ...
        @(x) ismember(x, {'proliferation', 'block_manual'}));
    addParameter(p, 'A1_matrix_unify_size', 'max_per_kernel', ...
        @(x) ismember(x, {'max_per_kernel', 'ref_slice'}));
    addParameter(p, 'window_type_proliferation', {'gaussian', 2.5});
    addParameter(p, 'interactive_size_adjust', false, @islogical);
    addParameter(p, 'use_matrix_format', true, @islogical);
    addParameter(p, 'show_reference_kernels', false, @islogical);
    parse(p, varargin{:});
    method = p.Results.block_init_method;
    A1_matrix_unify_size = p.Results.A1_matrix_unify_size;
    window_type_proliferation = p.Results.window_type_proliferation;
    interactive_size_adjust = p.Results.interactive_size_adjust;
    use_matrix_format = p.Results.use_matrix_format;
    show_reference_kernels = p.Results.show_reference_kernels;

    % Extract hierarchical → flat for reading
    if isfield(params, 'synGen')
        params = organizeParams(params, 'extract');
    end
    if isfield(data, 'synGen') || isfield(data, 'slice') || isfield(data, 'mcsbd_slice')
        data = organizeData(data, 'extract');
    end

    % Override method from params if set
    if isfield(params, 'block_init_method') && ~isempty(params.block_init_method)
        method = params.block_init_method;
    end

    % Store block params for write (so they appear in params.block)
    params.block_init_method = method;
    params.A1_matrix_unify_size = A1_matrix_unify_size;
    params.window_type_proliferation = window_type_proliferation;
    params.interactive_size_adjust = interactive_size_adjust;
    params.use_matrix_format = use_matrix_format;

    if strcmp(method, 'proliferation')
        % Delegate to proliferation wrapper (writes block init fields and calls organizeData/organizeParams)
        [data, params] = initializeProliferation(log, data, params, ...
            'window_type_proliferation', window_type_proliferation, ...
            'interactive_size_adjust', interactive_size_adjust, ...
            'use_matrix_format', use_matrix_format, ...
            'show_reference_kernels', show_reference_kernels, ...
            'A1_matrix_unify_size', A1_matrix_unify_size);
        % params already has block_init_method set; organizeParams was called inside initializeProliferation
        return;
    end

    % --- block_manual: no previous slice data needed; run manual selection, then propagate to all slices ---
    if ~strcmp(method, 'block_manual')
        error('block_init_method must be ''proliferation'' or ''block_manual''.');
    end

    % Require only what initializeKernelsRef does not validate: Y (all slices) and num_slices
    if ~isfield(data, 'Y')
        error('data must contain Y (block_manual uses all slices).');
    end
    if ~isfield(params, 'num_slices')
        error('params must contain num_slices (block_manual).');
    end
    num_kernels = params.num_kernels;
    num_slices  = params.num_slices;
    ref_slice   = params.ref_slice;

    % Always run manual selection for block_manual (avoids verifying params.slice.init_method)
    fprintf('Block init (block_manual): running manual ref-slice selection...\n');
    [data, params] = initializeKernelsRef(log, data, params);
    data = organizeData(data, 'extract');
    params = organizeParams(params, 'extract');
    A_ref = data.A_init;  % {1×K}

    % Kernel centers from ref: find position of each A_ref{k} in Y_ref via correlation
    kernel_centers = kernelCentersFromRef(data.Y_ref, A_ref, num_kernels);

    % Target kernel sizes: per-slice so each slice gets correct size at user-chosen centers (max_per_kernel)
    if isfield(params, 'target_kernel_sizes')
        use_slice_specific = isfield(params, 'target_kernel_size_type') && ...
            strcmp(params.target_kernel_size_type, 'kernel_sizes_all') && ndims(params.target_kernel_sizes) == 3;
    else
        use_slice_specific = false;
        if ~isfield(params, 'kernel_sizes')
            error('params must contain kernel_sizes or target_kernel_sizes (block_manual).');
        end
    end

    % Initialize kernels for all slices: same centers, per-slice kernel sizes (target_sizes_slice below)
    A1_all = cell(num_slices, num_kernels);
    fprintf('Initializing kernels for all %d slices (block_manual)...\n', num_slices);

    for s = 1:num_slices
        fprintf('  Slice %d/%d...', s, num_slices);
        if use_slice_specific
            target_sizes_slice = squeeze(params.target_kernel_sizes(s,:,:));
        elseif isfield(params, 'target_kernel_sizes')
            target_sizes_slice = params.target_kernel_sizes;
        else
            target_sizes_slice = reshape(params.kernel_sizes(s,:,:), [num_kernels, 2]);
        end
        A1_all(s,:) = initialize_kernels_proliferation(data.Y(:,:,s), num_kernels, ...
            kernel_centers, window_type_proliferation, target_sizes_slice, ...
            'interactive', interactive_size_adjust);
        fprintf(' done.\n');
    end

    % Visualize ref-slice kernels
    if show_reference_kernels
        figure('Name', 'IB01A: Block Init (block_manual) - Reference Slice');
        for k = 1:num_kernels
            subplot(1, num_kernels, k);
            imagesc(A1_all{ref_slice, k});
            colormap(gray);
            colorbar;
            title(sprintf('Kernel %d', k));
            axis square;
        end
        sgtitle(sprintf('Block Manual Init - Slice %d', ref_slice));
    end

    % Build A1_all_matrix with default A1_matrix_unify_size = 'max_per_kernel'
    if use_matrix_format
        A1_all_matrix = cell(num_kernels, 1);
        for k = 1:num_kernels
            switch A1_matrix_unify_size
                case 'max_per_kernel'
                    Hunif = max(arrayfun(@(s) size(A1_all{s,k}, 1), 1:num_slices));
                    Wunif = max(arrayfun(@(s) size(A1_all{s,k}, 2), 1:num_slices));
                case 'ref_slice'
                    Hunif = size(A1_all{ref_slice,k}, 1);
                    Wunif = size(A1_all{ref_slice,k}, 2);
                otherwise
                    Hunif = max(arrayfun(@(s) size(A1_all{s,k}, 1), 1:num_slices));
                    Wunif = max(arrayfun(@(s) size(A1_all{s,k}, 2), 1:num_slices));
            end
            A1_all_matrix{k} = zeros(Hunif, Wunif, num_slices);  % zero-padded; only kernel region filled per slice
            for s = 1:num_slices
                Ks = A1_all{s,k};
                [h, w] = size(Ks);
                r1_dest = floor((Hunif - h) / 2) + 1;
                c1_dest = floor((Wunif - w) / 2) + 1;
                r_dest = max(1, r1_dest) : min(Hunif, r1_dest + h - 1);
                c_dest = max(1, c1_dest) : min(Wunif, c1_dest + w - 1);
                r_src = (r_dest - r1_dest) + 1;
                c_src = (c_dest - c1_dest) + 1;
                A1_all_matrix{k}(r_dest, c_dest, s) = Ks(r_src, c_src);
            end
        end
        fprintf('Converted to matrix format: %d kernels, unify: %s\n\n', num_kernels, A1_matrix_unify_size);
    else
        A1_all_matrix = [];
    end

    fprintf('Block initialization (block_manual) complete.\n\n');

    % Write block init to flat data (organizeData maps to data.block)
    data.A1_all = A1_all;
    data.A1_all_matrix = A1_all_matrix;
    data.proliferation_kernel_centers = kernel_centers;

    LOGcomment = sprintf("Block init (block_manual): %d kernels, %d slices, centers from ref", ...
        num_kernels, num_slices);
    logUsedBlocks(log.path, log.file, "IB01A", LOGcomment, 0);

    data = organizeData(data, 'write');
    params = organizeParams(params, 'write');
end

%% Helper: kernel centers [K×2] from ref image and ref kernels (normxcorr2)
function centers = kernelCentersFromRef(Y_ref, A_ref, num_kernels)
    centers = zeros(num_kernels, 2);
    for k = 1:num_kernels
        template = A_ref{k};
        if isempty(template)
            error('A_ref{%d} is empty (block_manual).', k);
        end
        C = normxcorr2(template, Y_ref);
        [~, imax] = max(C(:));
        [pr, pc] = ind2sub(size(C), imax);
        % Peak in C corresponds to template center in image (normxcorr2 convention)
        ht = size(template, 1);
        wt = size(template, 2);
        center_r = pr - (ht - 1) / 2;
        center_c = pc - (wt - 1) / 2;
        centers(k, 1) = center_r;
        centers(k, 2) = center_c;
    end
end
