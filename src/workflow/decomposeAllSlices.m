function [data, params] = decomposeAllSlices(log, data, params, varargin)
%DECOMPOSEALLSLICES Wrapper for MT-SBD on all slices (DA01A).
%
%   Multi-kernel Tensor Shifted Blind Deconvolution across all energy slices,
%   using block initialization (IB01A) and optionally reference-slice X/b as init.
%
%   [data, params] = decomposeAllSlices(log, data, params, ...)
%
%   INPUTS:
%       log                 - Log struct with .path and .file
%       data                - Data from previous blocks (hierarchical or flat)
%       params              - Params from previous blocks (hierarchical or flat)
%
%       OPTIONAL (Name-Value), override params.block (full preset; [] = use params after extract):
%       Block-specific: 'use_reference_init', 'show_allslice_progress', 'maxIT_allslice'
%       Phase I: 'initial_iteration', 'lambda1'
%       Phase II: 'phase2_enable', 'lambda2', 'nrefine', 'kplus_factor'
%       Algorithm: 'signflip_threshold', 'xpos', 'getbias', 'Xsolve_method'
%       Init: 'use_xinit'. After extract, params.block.* overwrites params.slice.*; name-value overrides both.
%
%   OUTPUTS:
%       data  - Updated with data.block.Aout_all, Xout_all, bout_all, extras_all
%       params - Updated with params.block (full preset used for the run)
%
%   REQUIRED (before run): IB01A (data.block.A_all_init_matrix), DS01A (data.slice for ref init),
%       data.Y, data.X0; params.kernel_sizes, params.num_kernels, params.num_slices.
%
%   See also: MTSBD_synthetic_all_slice, decomposeReferenceSlice, organizeData, organizeParams

    % Validate inputs
    if ~isstruct(log) || ~isfield(log, 'path') || ~isfield(log, 'file')
        error('log must be a struct with .path and .file fields');
    end
    if ~isstruct(data) || ~isstruct(params)
        error('data and params must be structs');
    end

    % Parse optional name-value (full block preset; [] = use params after extract)
    p = inputParser;
    addParameter(p, 'use_reference_init', true, @islogical);
    addParameter(p, 'show_allslice_progress', true, @islogical);
    addParameter(p, 'maxIT_allslice', 15, @(x) isnumeric(x) && isscalar(x) && x > 0);
    addParameter(p, 'initial_iteration', []);
    addParameter(p, 'lambda1', []);
    addParameter(p, 'phase2_enable', []);
    addParameter(p, 'lambda2', []);
    addParameter(p, 'nrefine', []);
    addParameter(p, 'kplus_factor', []);
    addParameter(p, 'signflip_threshold', []);
    addParameter(p, 'xpos', []);
    addParameter(p, 'getbias', []);
    addParameter(p, 'Xsolve_method', '');
    addParameter(p, 'use_xinit', []);
    parse(p, varargin{:});
    nv = p.Results;

    % Extract hierarchical → flat (params.block overwrites params.slice for overlapping fields)
    if isfield(params, 'synGen')
        params = organizeParams(params, 'extract');
    end
    if isfield(data, 'synGen') || isfield(data, 'slice') || isfield(data, 'mcsbd_slice')
        data = organizeData(data, 'extract');
    end

    % Resolve block preset: name-value override else params (from block/slice)
    use_reference_init = nv.use_reference_init;
    if isfield(params, 'use_reference_init') && ~isempty(params.use_reference_init)
        use_reference_init = params.use_reference_init;
    end
    show_allslice_progress = nv.show_allslice_progress;
    if isfield(params, 'show_allslice_progress') && ~isempty(params.show_allslice_progress)
        show_allslice_progress = params.show_allslice_progress;
    end
    maxIT_allslice = round(nv.maxIT_allslice);
    if isfield(params, 'maxIT_allslice') && ~isempty(params.maxIT_allslice)
        maxIT_allslice = round(params.maxIT_allslice);
    end
    initial_iteration = params.initial_iteration;
    if ~isempty(nv.initial_iteration)
        initial_iteration = nv.initial_iteration;
    end
    lambda1 = params.lambda1;
    if ~isempty(nv.lambda1)
        lambda1 = nv.lambda1;
    end
    phase2_enable = params.phase2_enable;
    if ~isempty(nv.phase2_enable)
        phase2_enable = nv.phase2_enable;
    end
    lambda2 = params.lambda2;
    if ~isempty(nv.lambda2)
        lambda2 = nv.lambda2;
    end
    nrefine = params.nrefine;
    if ~isempty(nv.nrefine)
        nrefine = nv.nrefine;
    end
    kplus_factor = params.kplus_factor;
    if ~isempty(nv.kplus_factor)
        kplus_factor = nv.kplus_factor;
    end
    signflip_threshold = params.signflip_threshold;
    if ~isempty(nv.signflip_threshold)
        signflip_threshold = nv.signflip_threshold;
    end
    xpos = params.xpos;
    if ~isempty(nv.xpos)
        xpos = nv.xpos;
    end
    getbias = params.getbias;
    if ~isempty(nv.getbias)
        getbias = nv.getbias;
    end
    Xsolve_method = params.Xsolve_method;
    if ~isempty(nv.Xsolve_method)
        Xsolve_method = nv.Xsolve_method;
    end
    use_xinit = params.use_xinit;
    if ~isempty(nv.use_xinit)
        use_xinit = nv.use_xinit;
    end

    % Store block params for write (full preset used for the run)
    params.use_reference_init = use_reference_init;
    params.show_allslice_progress = show_allslice_progress;
    params.maxIT_allslice = maxIT_allslice;
    params.initial_iteration = initial_iteration;
    params.lambda1 = lambda1;
    params.phase2_enable = phase2_enable;
    params.lambda2 = lambda2;
    params.nrefine = nrefine;
    params.kplus_factor = kplus_factor;
    params.signflip_threshold = signflip_threshold;
    params.xpos = xpos;
    params.getbias = getbias;
    params.Xsolve_method = Xsolve_method;
    params.use_xinit = use_xinit;

    % Required data
    if ~isfield(data, 'Y')
        error('data must contain Y (run generation first).');
    end
    if ~isfield(data, 'A1_all_matrix')
        error('data must contain A1_all_matrix (run IB01A block init first).');
    end
    if ~isfield(data, 'X0')
        error('data must contain X0 (synGen).');
    end
    required_params = {'num_kernels', 'num_slices', 'kernel_sizes', 'initial_iteration', 'maxIT', ...
        'lambda1', 'phase2_enable', 'lambda2', 'nrefine', 'kplus_factor', 'signflip_threshold', ...
        'xpos', 'getbias', 'Xsolve_method'};
    for i = 1:length(required_params)
        if ~isfield(params, required_params{i})
            error('params must contain %s (run DS01A or set slice params).', required_params{i});
        end
    end

    num_kernels = params.num_kernels;
    num_slices  = params.num_slices;
    Y_used      = data.Y;
    A1_used     = data.A1_all_matrix;

    % lambda1, lambda2: core requires vector of length num_kernels; scalar → repmat, else error
    if isscalar(lambda1)
        lambda1 = repmat(lambda1, 1, num_kernels);
    elseif length(lambda1) ~= num_kernels
        error('lambda1 must be scalar or vector of length num_kernels (%d).', num_kernels);
    end
    if isscalar(lambda2)
        lambda2 = repmat(lambda2, 1, num_kernels);
    elseif length(lambda2) ~= num_kernels
        error('lambda2 must be scalar or vector of length num_kernels (%d).', num_kernels);
    end

    % Kernel sizes: max over slices per kernel [K×2]
    kernel_sizes_used = squeeze(max(params.kernel_sizes, [], 1));
    if size(kernel_sizes_used, 1) ~= num_kernels
        kernel_sizes_used = reshape(kernel_sizes_used, [num_kernels, 2]);
    end

    % Display functions
    if show_allslice_progress
        figure('Name', 'DA01A: All-Slice Decomposition Progress');
        dispfun_allslice = cell(1, num_kernels);
        X0_ref = data.X0_ref;
        for n = 1:num_kernels
            dispfun_allslice{n} = @(Y, A, X, kernel_sizes, kplus) ...
                showims(Y_used, A1_used{n}, X0_ref(:,:,n), A, X, kernel_sizes, kplus, 1);
        end
    else
        dispfun_allslice = cell(1, num_kernels);
        for n = 1:num_kernels
            dispfun_allslice{n} = @(Y, A, X, kernel_sizes, kplus) 0;
        end
    end

    % Prepare A0_used in matrix form (pad/crop to kernel_sizes_used, add noise per SNR)
    if isfield(data, 'synGen') && isfield(data.synGen, 'A0_noiseless')
        A0_noiseless = data.synGen.A0_noiseless;
    else
        A0_noiseless = data.A0_noiseless;
    end
    if isfield(params, 'synGen') && isfield(params.synGen, 'SNR')
        SNR_used = params.synGen.SNR;
    else
        SNR_used = params.SNR;
    end
    A0_used = padKernels(A0_noiseless, SNR_used, kernel_sizes_used);

    % Params for core: use resolved block preset (name-value or params.block/slice)
    params_allslice = struct();
    params_allslice.lambda1 = lambda1;
    params_allslice.phase2 = phase2_enable;
    params_allslice.kplus = ceil(kplus_factor * kernel_sizes_used);
    params_allslice.lambda2 = lambda2;
    params_allslice.nrefine = nrefine;
    params_allslice.signflip = signflip_threshold;
    params_allslice.xpos = xpos;
    params_allslice.getbias = getbias;
    params_allslice.Xsolve = Xsolve_method;
    % X0 from data.synGen, A0 = A0_used (matrix form, aligned with A1_all_matrix)
    if isfield(data, 'synGen') && isfield(data.synGen, 'X0')
        params_allslice.X0 = data.synGen.X0;
    else
        params_allslice.X0 = data.X0;
    end
    params_allslice.A0 = A0_used;

    % Optional: initialize X from reference slice (parse from data.slice or flat)
    if use_reference_init
        if isfield(data, 'slice') && isfield(data.slice, 'X') && isfield(data.slice, 'extras')
            X_ref = data.slice.X;
            extras_ref = data.slice.extras;
        elseif isfield(data, 'X') && isfield(data, 'extras')
            X_ref = data.X;
            extras_ref = data.extras;
        else
            error('use_reference_init requires data.slice (run DS01A first).');
        end
        params_allslice.xinit = cell(1, num_kernels);
        for k = 1:num_kernels
            params_allslice.xinit{k}.X = X_ref(:,:,k);
            b_temp = extras_ref.phase1.biter(k);
            params_allslice.xinit{k}.b = repmat(b_temp, [num_slices, 1]);
        end
    else
        params_allslice.xinit = [];
    end

    % LOG
    LOGcomment = sprintf("DA01A: use_ref_init=%d, maxIT=%d", use_reference_init, maxIT_allslice);
    logUsedBlocks(log.path, log.file, "DA01A", LOGcomment, 0);

    fprintf('Running MT-SBD decomposition on all %d slices...\n', num_slices);

    % Run all-slice decomposition
    tic;
    [Aout_slice, Xout_slice, bout_slice, slice_extras] = MTSBD_synthetic_all_slice(...
        Y_used, kernel_sizes_used, params_allslice, dispfun_allslice, ...
        A1_used, initial_iteration, maxIT_allslice);
    allslice_time = toc;

    fprintf('All-slice decomposition completed in %.2f seconds.\n', allslice_time);

    % Store block results in flat data (organizeData maps to data.block)
    data.Aout_all = Aout_slice;
    data.Xout_all = Xout_slice;
    data.bout_all = bout_slice;
    data.extras_all = slice_extras;

    LOGcomment = sprintf("Completed in %.2fs for %d slices", allslice_time, num_slices);
    logUsedBlocks(log.path, log.file, "  ^  ", LOGcomment, 0);

    data = organizeData(data, 'write');
    params = organizeParams(params, 'write');

    fprintf('Multi-slice decomposition complete.\n\n');
end
