function [data, params] = decomposeReferenceSlice(log, data, params, varargin)
%DECOMPOSEREFERENCESLICE Wrapper for MT-SBD on reference slice
%
%   Multi-kernel Tensor Shifted Blind Deconvolution for reference slice
%
%   [data, params, results] = decomposeReferenceSlice(log, data, params, ...)
%
%   INPUTS:
%       log                 - Log struct with .path and .file fields
%       data                - Data struct from previous blocks
%       params              - Parameter struct from previous blocks (hierarchical or flat)
%
%       OPTIONAL (Name-Value pairs):
%       Phase I settings:
%       'initial_iteration' - Manopt/FISTA inner iterations at start (default: 1)
%       'maxIT'             - Number of outer alternating iterations (default: 15)
%       'lambda1'           - L1 regularization for Phase I (default: auto-sized vector)
%
%       Phase II settings (refinement):
%       'phase2_enable'     - Enable Phase II refinement (default: false)
%       'lambda2'           - Final L1 regularization for Phase II (default: auto-sized)
%       'nrefine'           - Number of refinement steps (default: 5)
%       'kplus_factor'      - Sphere lifting padding factor (default: 0.5)
%
%       Algorithm parameters:
%       'signflip_threshold' - Sign flip detection threshold (default: 0.2)
%       'xpos'              - Enforce positive activations (default: true)
%       'getbias'           - Extract constant bias term (default: true)
%       'Xsolve_method'     - Solver: 'FISTA' or 'pdNCG' (default: 'FISTA')
%
%       Initialization options:
%       'use_xinit'         - Initial X guess, [] for none (default: [])
%
%       Display options:
%       'show_progress'     - Show MT-SBD optimization progress (default: true)
%
%   OUTPUTS:
%       data                - Updated data struct with new fields in mcsbd_slice:
%                             data.mcsbd_slice.A - Deconvolved kernels
%                             data.mcsbd_slice.X - Activation maps
%                             data.mcsbd_slice.b - Bias terms
%                             data.mcsbd_slice.extras - Optimization info
%                             data.mcsbd_slice.mtsbd_time - Execution time
%                             data.mcsbd_slice.final_metrics - Activation metrics
%                             data.mcsbd_slice.final_kernel_quality - Kernel quality
%       params              - Updated parameter struct with MT-SBD settings
%       results             - Results struct (for backward compatibility, contains same as data.mcsbd_slice)
%
%   DESCRIPTION:
%       This wrapper encapsulates the MT-SBD (Multi-kernel Tensor Shifted
%       Blind Deconvolution) process for the reference slice. It handles:
%       - Setting up display functions for monitoring
%       - Packaging parameters for MTSBD_synthetic
%       - Running the MT-SBD algorithm
%       - Computing and displaying quality metrics
%       - Storing results in organized structures
%
%   EXAMPLE:
%       % Standard MT-SBD
%       [data, params, results] = decomposeReferenceSlice(log, data, params);
%
%       % With Phase II refinement
%       [data, params, results] = decomposeReferenceSlice(log, data, params, ...
%           'phase2_enable', true, 'maxIT', 20);
%
%   See also: MTSBD_synthetic, Asolve_Manopt_tunable, Xsolve_FISTA_tunable

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
    
    % Extract parameters from hierarchical structure (if needed)
    if isfield(params, 'synGen')
        params = organizeParams(params, 'extract');
    end
    
    % Extract data from hierarchical structure (if needed)
    if isfield(data, 'synGen') || isfield(data, 'mcsbd_slice')
        data = organizeData(data, 'extract');
    end
    
    % All parameters are required; do not give default values
    required_fields = {'initial_iteration', 'maxIT', 'lambda1', 'phase2_enable', ...
        'lambda2', 'nrefine', 'kplus_factor', 'signflip_threshold', 'xpos', ...
        'getbias', 'Xsolve_method', 'use_xinit', 'show_progress'};
    for i = 1:length(required_fields)
        if ~isfield(params, required_fields{i})
            error('Parameter "%s" is required in params struct.', required_fields{i});
        end
    end

    initial_iteration   = params.initial_iteration;
    maxIT               = params.maxIT;
    lambda1             = params.lambda1;
    phase2_enable       = params.phase2_enable;
    lambda2             = params.lambda2;
    nrefine             = params.nrefine;
    kplus_factor        = params.kplus_factor;
    signflip_threshold  = params.signflip_threshold;
    xpos                = params.xpos;
    getbias             = params.getbias;
    Xsolve_method       = params.Xsolve_method;
    use_xinit           = params.use_xinit;
    show_progress       = params.show_progress;
    
    % Validate Xsolve_method
    if ~ismember(Xsolve_method, {'FISTA', 'pdNCG'})
        error('Xsolve_method must be ''FISTA'' or ''pdNCG''');
    end
    
    % Validate required fields in data struct (flat structure, after unpacking)
    if ~isfield(data, 'Y_ref')
        error('Data struct must contain field: Y_ref');
    end
    if ~isfield(data, 'X0_ref')
        error('Data struct must contain field: X0_ref');
    end
    if ~isfield(data, 'A0_ref')
        error('Data struct must contain field: A0_ref');
    end
    if ~isfield(data, 'A_init')
        error('Data struct must contain field: A_init (from initializeKernelsRef or autoInitializeKernels)');
    end
    
    % Validate required fields in params struct
    required_params = {'num_kernels', 'ref_slice', 'kernel_sizes'};
    for i = 1:length(required_params)
        if ~isfield(params, required_params{i})
            error('Params struct must contain field: %s', required_params{i});
        end
    end
    
    % Auto-size lambda vectors if not provided
    if isempty(lambda1)
        lambda1 = repmat(3e-2, 1, params.num_kernels);
    elseif isscalar(lambda1)
        lambda1 = repmat(lambda1, 1, params.num_kernels);
    elseif length(lambda1) ~= params.num_kernels
        error('lambda1 must be scalar or vector of length num_kernels (%d)', params.num_kernels);
    end
    
    if isempty(lambda2)
        lambda2 = repmat(1e-2, 1, params.num_kernels);
    elseif isscalar(lambda2)
        lambda2 = repmat(lambda2, 1, params.num_kernels);
    elseif length(lambda2) ~= params.num_kernels
        error('lambda2 must be scalar or vector of length num_kernels (%d)', params.num_kernels);
    end
    
    % LOG: function start
    LOGcomment = sprintf("DS01A: maxIT=%d, lambda1=%s, phase2=%d", maxIT, mat2str(lambda1, 3), phase2_enable);
    LOGcomment = logUsedBlocks(log.path, log.file, "DS01A", LOGcomment, 0);
    
    % Set up display functions for monitoring
    fprintf('  Setting up MT-SBD...\n');
    if show_progress
        figure('Name', 'MT-SBD Progress');
        dispfun = cell(1, params.num_kernels);
        for n = 1:params.num_kernels
            dispfun{n} = @(Y, A, X, kernel_sizes, kplus) ...
                showims(data.Y_ref, data.A_init{n}, data.X0_ref(:,:,n), A, X, kernel_sizes, kplus, 1);
        end
    else
        dispfun = cell(1, params.num_kernels);
        for n = 1:params.num_kernels
            dispfun{n} = @(Y, A, X, kernel_sizes, kplus) 0;
        end
    end
    
    % Package parameters
    kernel_sizes_ref = reshape(params.kernel_sizes(params.ref_slice,:,:), [params.num_kernels, 2]);
    sbd_params = struct();
    sbd_params.lambda1 = lambda1;
    sbd_params.phase2 = phase2_enable;
    sbd_params.kplus = ceil(kplus_factor * kernel_sizes_ref);
    sbd_params.lambda2 = lambda2;
    sbd_params.nrefine = nrefine;
    sbd_params.signflip = signflip_threshold;
    sbd_params.xpos = xpos;
    sbd_params.getbias = getbias;
    sbd_params.Xsolve = Xsolve_method;
    sbd_params.X0 = data.X0_ref;
    sbd_params.A0 = data.A0_ref;
    sbd_params.xinit = use_xinit;
    
    % Run MT-SBD on reference slice
    tic;
    [A_ref, X_ref, bout, extras] = MTSBD_synthetic(...
        data.Y_ref, kernel_sizes_ref, sbd_params, dispfun, data.A_init, initial_iteration, maxIT);
    mtsbd_time = toc;
    
    fprintf('  MT-SBD completed in %.2f seconds.\n', mtsbd_time);
    
    % Display final quality metrics
    fprintf('\n  Final Quality Metrics:\n');
    final_metrics = extras.phase1.activation_metrics(end,:);
    final_kernel_quality = extras.phase1.kernel_quality_factors(end,:);
    for k = 1:params.num_kernels
        fprintf('    Kernel %d - Activation: %.4f, Quality: %.4f\n', ...
            k, final_metrics(k), final_kernel_quality(k));
    end
    fprintf('\n');
    
    % LOG: MT-SBD results
    LOGcomment = sprintf("Completed in %.2fs, Final metrics: Act=%s, Qual=%s", ...
        mtsbd_time, mat2str(final_metrics, 3), mat2str(final_kernel_quality, 3));
    LOGcomment = logUsedBlocks(log.path, log.file, "  ^  ", LOGcomment, 0);
    
    % Store results in data
    data.X = X_ref;
    data.b = bout;
    data.extras = extras;
    data.mtsbd_time = mtsbd_time;
    data.final_metrics = final_metrics;
    data.final_kernel_quality = final_kernel_quality;
    
    % Update params struct with MT-SBD settings (flat structure)
    params.initial_iteration = initial_iteration;
    params.maxIT = maxIT;
    params.lambda1 = lambda1;
    params.phase2_enable = phase2_enable;
    params.lambda2 = lambda2;
    params.nrefine = nrefine;
    params.kplus_factor = kplus_factor;
    params.signflip_threshold = signflip_threshold;
    params.xpos = xpos;
    params.getbias = getbias;
    params.Xsolve_method = Xsolve_method;
    params.kernel_sizes_ref = kernel_sizes_ref;
    
    % Write data and params to hierarchical structure for storage
    data = organizeData(data, 'write');
    params = organizeParams(params, 'write');
    
    fprintf('  Reference slice MT-SBD complete.\n');
end


