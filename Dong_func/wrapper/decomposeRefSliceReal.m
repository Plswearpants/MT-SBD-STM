function [log, data, params, meta, cfg] = decomposeRefSliceReal(log, data, params, meta, cfg)
%DECOMPOSEREFSLICEREAL Decompose reference slice for real data (Block 4).
%
%   [log, data, params, meta, cfg] = decomposeRefSliceReal(log, data, params, meta, cfg)
%
%   This wrapper encapsulates the "Pick reference slice" + "Block 4: Find
%   Optimal Activation for Reference Slice" logic from
%   scripts/MTSBD_block_realdata1.m. It:
%       - lets the user choose a reference slice and number of kernels
%       - initializes kernels on the reference slice (initialize_kernels)
%       - enforces kernel polarity
%       - estimates noise
%       - runs MT_SBD on the reference slice
%       - stores results under data.real and params.refSlice
%
%   Retired pieces (SNR computation, padded reference run, manual X
%   comparison) are intentionally omitted.

    arguments
        log  struct
        data struct
        params struct
        meta struct
        cfg  struct
    end

    if ~isfield(data, "real") || ~isfield(data.real, "Y")
        error('decomposeRefSliceReal: data.real.Y is missing. Run preprocessRealData first.');
    end

    Y = data.real.Y;
    rangetype = 'dynamic';

    % ---------------------------------------------------------------------
    % 4_pre: Pick reference slice and number of kernels
    % ---------------------------------------------------------------------
    figure;
    d3gridDisplay(Y, rangetype);
    params.refSlice.ref_slice = input('Enter reference slice number: ');
    num_kernels = input('enter the number of kernels you wish to apply: ');

    if isempty(params.refSlice.ref_slice) && isfield(cfg, "reference") && ~isempty(cfg.reference.default_ref_slice)
        params.refSlice.ref_slice = cfg.reference.default_ref_slice;
    end
    if isempty(num_kernels) && isfield(cfg, "reference") && ~isempty(cfg.reference.default_num_kernels)
        num_kernels = cfg.reference.default_num_kernels;
    end
    close;

    % Validate input
    if isempty(params.refSlice.ref_slice) || ~isnumeric(params.refSlice.ref_slice) || ...
            params.refSlice.ref_slice < 1 || params.refSlice.ref_slice > size(Y, 3)
        error('Invalid reference slice. Please enter a number between 1 and %d.', size(Y, 3));
    end
    fprintf('Using slice %d as reference\n', params.refSlice.ref_slice);

    % Normalize and project
    Y = normalizeBackgroundToZeroMean3D(Y, 'dynamic', params.refSlice.ref_slice);
    Y = proj2oblique(Y);

    % Extract reference slice and normalize
    Y_ref = Y(:,:,params.refSlice.ref_slice);
    Y_ref = normalizeBackgroundToZeroMean3D(Y_ref, 'dynamic');
    Y_ref = proj2oblique(Y_ref);

    % Display reference slice
    figure;
    imagesc(Y_ref);
    colorbar;
    title(sprintf('Reference Slice %d', params.refSlice.ref_slice));
    axis square;
    close;

    % ---------------------------------------------------------------------
    % Initialize reference kernels (initialize_kernels)
    % ---------------------------------------------------------------------
    same_size   = cfg.reference.same_size;
    kerneltype  = cfg.reference.kerneltype;
    window_type = cfg.reference.window_type;

    if same_size
        square_size = cfg.reference.square_size;
        kernel_sizes = repmat(square_size, [num_kernels, 1]);
        [A1_ref, A1_ref_crop] = initialize_kernels(Y_ref, num_kernels, kernel_sizes, kerneltype, window_type);
    else
        A1_ref = cell(1, num_kernels);
        A1_ref_crop = cell(1, num_kernels);
        kernel_sizes = zeros(num_kernels, 2);
        for k = 1:num_kernels
            fprintf('Select region for kernel %d/%d\n', k, num_kernels);
            [square_size, position, mask] = squareDrawSize(Y_ref); %#ok<NASGU>
            [A1_ref{k}, ~] = gridCropMask(Y_ref, mask);
            A1_ref_crop{k} = A1_ref{k};
            A1_ref{k} = proj2oblique(A1_ref{k});
            kernel_sizes(k,:) = size(A1_ref_crop{k});
        end
    end

    % Enforce kernel polarity
    for k = 1:num_kernels
        [A1_ref{k}, flipped_ref] = enforce_kernel_polarity(A1_ref{k}, A1_ref_crop{k});
        if flipped_ref
            fprintf('[kernel flip] reference kernel %d flipped\n', k);
        end
    end

    % Quick visualization of initialized kernels
    figure;
    for k = 1:num_kernels
        subplot(1, num_kernels, k);
        imagesc(A1_ref{k});
        colorbar;
        axis square;
    end

    % ---------------------------------------------------------------------
    % Noise estimation and MT_SBD on reference slice
    % ---------------------------------------------------------------------
    eta_data = estimate_noise(Y_ref, 'std');

    figure;
    dispfun = cell(1, num_kernels);
    for n = 1:num_kernels
        dispfun{n} = @(Y_, A, X, kernel_sizes_, kplus) ... %#ok<NASGU,INUSD>
            showims(Y_ref, A1_ref{n}, X, A, X, kernel_sizes_, kplus, 1);
    end

    miniloop_iteration = cfg.sliceRun.miniloop_iteration;
    outerloop_maxIT    = cfg.sliceRun.outerloop_maxIT;

    params_ref = struct();
    params_ref.lambda1  = cfg.sliceRun.lambda1;
    params_ref.phase2   = cfg.sliceRun.phase2;
    params_ref.kplus    = ceil(cfg.sliceRun.kplus_factor * kernel_sizes);
    params_ref.lambda2  = cfg.sliceRun.lambda2;
    params_ref.nrefine  = cfg.sliceRun.nrefine;
    params_ref.signflip = cfg.sliceRun.signflip;
    params_ref.xpos     = cfg.sliceRun.xpos;
    params_ref.getbias  = cfg.sliceRun.getbias;
    params_ref.Xsolve   = cfg.sliceRun.Xsolve;
    params_ref.noise_var = eta_data;

    [A_ref, X_ref, b_ref, extras_ref] = MT_SBD( ...
        Y_ref, kernel_sizes, params_ref, dispfun, A1_ref, ...
        miniloop_iteration, outerloop_maxIT);

    % Visualize reference result (legacy helper)
    [Y_rec_ref, Y_rec_all_ref] = visualizeRealResult(Y_ref, A_ref, X_ref, b_ref, extras_ref); %#ok<NASGU>

    % ---------------------------------------------------------------------
    % Store results
    % ---------------------------------------------------------------------
    if ~isfield(data.real, "ref"); data.real.ref = struct(); end

    data.real.ref.Y_ref     = Y_ref;
    data.real.ref.A1_ref    = A1_ref;
    data.real.ref.A1_ref_crop = A1_ref_crop;
    data.real.ref.kernel_sizes = kernel_sizes;
    data.real.ref.A_ref     = A_ref;
    data.real.ref.X_ref     = X_ref;
    data.real.ref.b_ref     = b_ref;
    data.real.ref.extras_ref = extras_ref;
    data.real.Y             = Y;      % normalized 3D volume

    params.refSlice.num_kernels = num_kernels;
    params.refSlice.kernel_sizes = kernel_sizes;
    params.refSlice.params_ref   = params_ref;

    meta.stage = "pre-run";

end

