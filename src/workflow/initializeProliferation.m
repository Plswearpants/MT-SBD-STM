function [data, params] = initializeProliferation(log, data, params, varargin)
%INITIALIZEPROLIFERATION Wrapper for initializing kernels across all slices (proliferation)
%
%   Uses the most isolated points from the reference slice as centers to
%   initialize kernels on every energy slice. This "proliferation" ensures
%   consistent spatial positioning across the energy dimension.
%
%   [data, params] = initializeProliferation(log, data, params, ...)
%
%   INPUTS:
%       log                 - Log struct with .path and .file fields
%       data                - Data struct from previous blocks (must contain
%                             data.slice.most_isolated_points from IS01A; or data.most_isolated_points after extract)
%       params              - Parameter struct from previous blocks (hierarchical or flat)
%
%       OPTIONAL (Name-Value pairs):
%       'window_type_proliferation' - Window function for proliferated kernels
%                                     (default: {'gaussian', 2.5})
%       'interactive_size_adjust'   - Allow interactive kernel size adjustment (default: false)
%       'use_matrix_format'        - Output A1_all_matrix as {K×1} 3D arrays (default: true)
%       'show_reference_kernels'   - Display initialized kernels for reference slice (default: true)
%       'A1_matrix_unify_size'     - How to unify sizes when building A1_all_matrix (sizes can
%                                     differ per slice). Only used when use_matrix_format is true.
%                                     (default: 'max_per_kernel')
%         'max_per_kernel' - Pad each slice to max height/width over all slices for that kernel.
%         'ref_slice'      - Use reference-slice dimensions; pad/crop other slices to match.
%
%   OUTPUTS:
%       data                - Updated data struct with new fields in block (organizeData maps to block):
%                             data.block.A_all_init - Cell {S×K} theoretical kernel per slice (sizes may differ)
%                             data.block.A_all_init_matrix - Cell {K×1} [H×W×S] unified size (if use_matrix_format);
%                               H,W chosen by A1_matrix_unify_size (default: max per kernel, center-padded)
%                             data.block.init_kernel_centers - [K×2] centers used (proliferation)
%       params              - Unchanged (no new params from this block)
%
%   REQUIRED (in data/params after extract):
%       data.Y, data.most_isolated_points
%       params.num_kernels, params.num_slices, params.ref_slice
%       params.kernel_sizes and/or params.target_kernel_sizes, params.target_kernel_size_type
%
%   See also: initialize_kernels_proliferation, findIsolatedPoints, formatWindowType

    % Validate required inputs
    if ~isstruct(log) || ~isfield(log, 'path') || ~isfield(log, 'file')
        error('log must be a struct with .path and .file fields');
    end
    if ~isstruct(data)
        error('data must be a struct');
    end
    if ~isstruct(params)
        error('params must be a struct');
    end

    % Parse optional name-value pairs
    p = inputParser;
    addParameter(p, 'window_type_proliferation', {'gaussian', 2.5});
    addParameter(p, 'interactive_size_adjust', false, @islogical);
    addParameter(p, 'use_matrix_format', true, @islogical);
    addParameter(p, 'show_reference_kernels', true, @islogical);
    addParameter(p, 'A1_matrix_unify_size', 'max_per_kernel', ...
        @(x) ismember(x, {'max_per_kernel', 'ref_slice'}));
    parse(p, varargin{:});
    window_type_proliferation = p.Results.window_type_proliferation;
    interactive_size_adjust  = p.Results.interactive_size_adjust;
    use_matrix_format       = p.Results.use_matrix_format;
    show_reference_kernels   = p.Results.show_reference_kernels;
    A1_matrix_unify_size    = p.Results.A1_matrix_unify_size;

    % Extract from hierarchical structure (if needed)
    if isfield(params, 'synGen')
        params = organizeParams(params, 'extract');
    end
    if isfield(data, 'synGen') || isfield(data, 'slice') || isfield(data, 'mcsbd_slice')
        data = organizeData(data, 'extract');
    end

    % Validate required fields
    if ~isfield(data, 'Y')
        error('data must contain Y (observation stack)');
    end
    if ~isfield(data, 'most_isolated_points')
        error('data must contain most_isolated_points (run IS01A first)');
    end
    for f = {'num_kernels', 'num_slices', 'ref_slice'}
        if ~isfield(params, f{1})
            error('params must contain %s', f{1});
        end
    end

    num_kernels = params.num_kernels;
    num_slices  = params.num_slices;
    ref_slice   = params.ref_slice;

    % Build kernel_centers matrix from most_isolated_points
    kernel_centers = zeros(num_kernels, 2);
    for k = 1:num_kernels
        if ~isempty(data.most_isolated_points{k})
            kernel_centers(k,:) = data.most_isolated_points{k};
        else
            error('No isolated point found for kernel %d', k);
        end
    end

    % Determine target kernel sizes per slice (from IS01A: params.target_kernel_sizes)
    if ~isfield(params, 'target_kernel_sizes')
        error('params must contain target_kernel_sizes (run IS01A first)');
    end
    use_slice_specific = isfield(params, 'target_kernel_size_type') && ...
        strcmp(params.target_kernel_size_type, 'kernel_sizes_all') && ndims(params.target_kernel_sizes) == 3;

    % Initialize kernels for all slices
    A1_all = cell(num_slices, num_kernels);
    fprintf('Initializing kernels for all %d slices using proliferation...\n', num_slices);

    for s = 1:num_slices
        fprintf('  Slice %d/%d...', s, num_slices);
        if use_slice_specific
            target_sizes_slice = squeeze(params.target_kernel_sizes(s,:,:));
        else
            target_sizes_slice = params.target_kernel_sizes;
        end
        A1_all(s,:) = initialize_kernels_proliferation(data.Y(:,:,s), num_kernels, ...
            kernel_centers, window_type_proliferation, target_sizes_slice, ...
            'interactive', interactive_size_adjust);
        fprintf(' done.\n');
    end

    % Visualize initialized kernels for reference slice
    if show_reference_kernels
        figure('Name', 'IP01A: Initialized Kernels (Reference Slice)');
        for k = 1:num_kernels
            subplot(1, num_kernels, k);
            imagesc(A1_all{ref_slice, k});
            colormap(gray);
            colorbar;
            title(sprintf('Initialized Kernel %d', k));
            axis square;
        end
        sgtitle(sprintf('Proliferated Kernels - Slice %d', ref_slice));
    end

    % Convert to matrix format if requested: unify sizes per A1_matrix_unify_size
    % A1_all keeps theoretical (slice-specific) sizes; A1_all_matrix uses unified [H x W x S].
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
            A1_all_matrix{k} = zeros(Hunif, Wunif, num_slices);  % zero-filled; only kernel region written below
            for s = 1:num_slices
                Ks = A1_all{s,k};
                [h, w] = size(Ks);
                % Center Ks in [Hunif x Wunif] frame: pad with zeros if smaller, crop if larger
                r1_dest = floor((Hunif - h) / 2) + 1;
                c1_dest = floor((Wunif - w) / 2) + 1;
                r_dest = max(1, r1_dest) : min(Hunif, r1_dest + h - 1);
                c_dest = max(1, c1_dest) : min(Wunif, c1_dest + w - 1);
                r_src = (r_dest - r1_dest) + 1;
                c_src = (c_dest - c1_dest) + 1;
                A1_all_matrix{k}(r_dest, c_dest, s) = Ks(r_src, c_src);
            end
        end
        fprintf('Converted to matrix format: %d kernels x [%dx%dx%d] (unify: %s)\n\n', ...
            num_kernels, size(A1_all_matrix{1},1), size(A1_all_matrix{1},2), num_slices, A1_matrix_unify_size);
    else
        A1_all_matrix = [];
    end

    fprintf('Kernel initialization complete for all slices.\n\n');

    % Store in flat data struct
    data.A1_all = A1_all;
    data.A1_all_matrix = A1_all_matrix;
    data.proliferation_kernel_centers = kernel_centers;

    % LOG: completion (use formatWindowType for window string)
    window_str = formatWindowType(window_type_proliferation);
    LOGcomment = sprintf("Initialized %d kernels for %d slices at centers: %s", ...
        num_kernels, num_slices, mat2str(kernel_centers));
    logUsedBlocks(log.path, log.file, "  ^  ", LOGcomment, 0);

    % Convert to hierarchical for storage
    data = organizeData(data, 'write');
    params = organizeParams(params, 'write');
end
