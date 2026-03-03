function [log, data, params, meta, cfg] = initProliferationReal(log, data, params, meta, cfg)
%INITPROLIFERATIONREAL Initialize kernels for all slices (proliferation).
%
%   [log, data, params, meta, cfg] = initProliferationReal(log, data, params, meta, cfg)
%
%   This wrapper encapsulates the "Block 3: Find Most Isolated Points and
%   Initialize Kernels" logic from scripts/MTSBD_block_realdata1.m, with
%   the most-isolated-points AUTO mode treated as retired. It:
%       - lets the user manually select kernel centers on Y_ref
%       - proliferates kernels across all slices via initialize_kernels_proliferation
%       - enforces kernel polarity per slice
%       - converts A1_all to matrix form A1_all_matrix
%       - estimates per-slice noise (eta_data3d)
%
%   Results are stored under data.real.proliferation.* and params.proliferation.

    arguments
        log  struct
        data struct
        params struct
        meta struct
        cfg  struct
    end

    if ~isfield(data, "real") || ~isfield(data.real, "Y")
        error('initProliferationReal: data.real.Y is missing. Run preprocessRealData first.');
    end
    if ~isfield(data.real, "ref") || ~isfield(data.real.ref, "Y_ref")
        error('initProliferationReal: data.real.ref.Y_ref is missing. Run decomposeRefSliceReal first.');
    end

    Y     = data.real.Y;
    Y_ref = data.real.ref.Y_ref;
    kernel_sizes = data.real.ref.kernel_sizes;
    num_kernels  = params.refSlice.num_kernels;
    ref_slice    = params.refSlice.ref_slice;

    num_slices = size(Y, 3);

    % ---------------------------------------------------------------------
    % Manual kernel center selection (AUTO isolation logic retired)
    % ---------------------------------------------------------------------
    fprintf('Manual kernel center selection mode...\n');

    figure('Name', 'Manual Kernel Center Selection (Proliferation)');
    imagesc(Y_ref);
    axis square;
    title('Click on centers to select kernel positions. Press Enter when done.');
    colormap(gray);
    colorbar;

    kernel_centers = zeros(num_kernels, 2);
    for k = 1:num_kernels
        fprintf('Click on center for kernel %d/%d\n', k, num_kernels);
        [x, y] = ginput(1);
        kernel_centers(k,:) = [round(y), round(x)];  % [row, col]

        hold on;
        scatter(x, y, 100, 'r', '*');
        text(x+5, y+5, sprintf('K%d', k), 'Color', 'red', ...
            'FontSize', 12, 'FontWeight', 'bold');
        hold off;
    end

    fprintf('Kernel centers selected:\n');
    for k = 1:num_kernels
        fprintf('Kernel %d: (%d, %d)\n', k, kernel_centers(k,1), kernel_centers(k,2));
    end

    % Target kernel sizes: use reference kernel sizes by default
    target_kernel_sizes = kernel_sizes;

    % ---------------------------------------------------------------------
    % Initialize kernels for all slices (initialize_kernels_proliferation)
    % ---------------------------------------------------------------------
    A1_all      = cell(num_slices, num_kernels);
    A1_all_crop = cell(num_slices, num_kernels);

    matrix      = cfg.blockInit.use_matrix;
    change_size = cfg.blockInit.change_size;
    window_type = cfg.reference.window_type;

    for s = 1:num_slices
        fprintf('Initializing kernels for slice %d/%d...\n', s, num_slices);
        if matrix
            [A1_all(s,:), A1_all_crop(s,:)] = initialize_kernels_proliferation( ...
                Y(:,:,s), num_kernels, kernel_centers, window_type, ...
                target_kernel_sizes, 'interactive', change_size);
        else
            [A1_all(s,:), A1_all_crop(s,:)] = initialize_kernels_proliferation( ...
                Y(:,:,s), num_kernels, kernel_centers, window_type, ...
                squeeze(kernel_sizes(s,:,:)), 'interactive', change_size);
        end
    end

    % Enforce kernel polarity per slice
    flip_slices_by_kernel = cell(1, num_kernels);
    for s = 1:num_slices
        for k = 1:num_kernels
            [A1_all{s,k}, flipped] = enforce_kernel_polarity(A1_all{s,k}, A1_all_crop{s,k});
            if flipped
                flip_slices_by_kernel{k}(end+1) = s; %#ok<AGROW>
                fprintf('[kernel flip] kernel %d flipped at slice %d\n', k, s);
            end
        end
    end

    fprintf('==== Kernel flip summary (block init) ====\n');
    for k = 1:num_kernels
        if isempty(flip_slices_by_kernel{k})
            fprintf('Kernel %d: flipped slices = []\n', k);
        else
            flip_slices_by_kernel{k} = unique(flip_slices_by_kernel{k});
            fprintf('Kernel %d: flipped slices = %s\n', k, mat2str(flip_slices_by_kernel{k}));
        end
    end

    % Visualize initialized kernels for reference slice
    figure('Name', 'Initialized Kernels (Reference Slice)');
    for k = 1:num_kernels
        subplot(1, num_kernels, k);
        imagesc(A1_all{ref_slice,k});
        colormap(gray);
        colorbar;
        title(sprintf('Initialized Kernel %d', k));
        axis square;
    end

    % Convert A1_all to matrix form
    A1_all_matrix = cell(num_kernels, 1);
    for k = 1:num_kernels
        A1_all_matrix{k} = zeros(size(A1_all{1,k},1), size(A1_all{1,k},2), num_slices);
        for s = 1:num_slices
            A1_all_matrix{k}(:,:,s) = A1_all{s,k};
        end
    end

    % Noise estimate per slice
    eta_data3d = estimate_noise3D(Y, 'std');

    % ---------------------------------------------------------------------
    % Store results
    % ---------------------------------------------------------------------
    if ~isfield(data.real, "proliferation")
        data.real.proliferation = struct();
    end

    data.real.proliferation.A1_all       = A1_all;
    data.real.proliferation.A1_all_crop  = A1_all_crop;
    data.real.proliferation.A1_all_matrix = A1_all_matrix;
    data.real.proliferation.kernel_centers = kernel_centers;
    data.real.proliferation.target_kernel_sizes = target_kernel_sizes;
    data.real.proliferation.eta_data3d   = eta_data3d;

    params.proliferation.kernel_centers      = kernel_centers;
    params.proliferation.target_kernel_sizes = target_kernel_sizes;

    % Stage remains "pre-run" (still preparation for block run)
    meta.stage = "pre-run";

end

