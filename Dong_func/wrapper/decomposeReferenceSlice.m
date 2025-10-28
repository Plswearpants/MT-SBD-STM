function [data, params, results] = decomposeReferenceSlice(data, params, varargin)
%DECOMPOSEREFERENCESLICE Wrapper for MT-SBD on reference slice
%
%   Multi-kernel Tensor Shifted Blind Deconvolution for reference slice
%
%   [data, params, results] = decomposeReferenceSlice(data, params, ...)
%
%   INPUTS:
%       data                - Data struct from previous blocks
%       params              - Parameter struct from previous blocks
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
%       data                - Updated data struct with new fields:
%                             data.A_ref  - Deconvolved kernels
%                             data.X_ref  - Activation maps
%                             data.b_ref  - Bias terms
%       params              - Updated parameter struct with MT-SBD settings
%       results             - Results struct containing:
%                             results.A - Refined kernels
%                             results.X - Activation maps
%                             results.b - Bias terms
%                             results.extras - Detailed optimization info
%                             results.params - Algorithm parameters used
%                             results.mtsbd_time - Execution time
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
%       [data, params, results] = decomposeReferenceSlice(data, params);
%
%       % With Phase II refinement
%       [data, params, results] = decomposeReferenceSlice(data, params, ...
%           'phase2_enable', true, 'maxIT', 20);
%
%   See also: MTSBD_synthetic, Asolve_Manopt_tunable, Xsolve_FISTA_tunable

    % Parse input arguments
    p = inputParser;
    addRequired(p, 'data', @isstruct);
    addRequired(p, 'params', @isstruct);
    
    % Phase I settings
    addParameter(p, 'initial_iteration', 1, @isnumeric);
    addParameter(p, 'maxIT', 15, @isnumeric);
    addParameter(p, 'lambda1', [], @(x) isempty(x) || isnumeric(x));
    
    % Phase II settings
    addParameter(p, 'phase2_enable', false, @islogical);
    addParameter(p, 'lambda2', [], @(x) isempty(x) || isnumeric(x));
    addParameter(p, 'nrefine', 5, @isnumeric);
    addParameter(p, 'kplus_factor', 0.5, @isnumeric);
    
    % Algorithm parameters
    addParameter(p, 'signflip_threshold', 0.2, @isnumeric);
    addParameter(p, 'xpos', true, @islogical);
    addParameter(p, 'getbias', true, @islogical);
    addParameter(p, 'Xsolve_method', 'FISTA', @(x) ismember(x, {'FISTA', 'pdNCG'}));
    
    % Initialization options
    addParameter(p, 'use_xinit', [], @(x) isempty(x) || isstruct(x));
    
    % Display options
    addParameter(p, 'show_progress', true, @islogical);
    
    parse(p, data, params, varargin{:});
    
    % Extract parameters
    initial_iteration = p.Results.initial_iteration;
    maxIT = p.Results.maxIT;
    lambda1 = p.Results.lambda1;
    phase2_enable = p.Results.phase2_enable;
    lambda2 = p.Results.lambda2;
    nrefine = p.Results.nrefine;
    kplus_factor = p.Results.kplus_factor;
    signflip_threshold = p.Results.signflip_threshold;
    xpos = p.Results.xpos;
    getbias = p.Results.getbias;
    Xsolve_method = p.Results.Xsolve_method;
    use_xinit = p.Results.use_xinit;
    show_progress = p.Results.show_progress;
    
    % Validate required fields in data struct
    required_fields = {'Y_ref', 'A1', 'X0_ref', 'A0_ref'};
    for i = 1:length(required_fields)
        if ~isfield(data, required_fields{i})
            error('Data struct must contain field: %s', required_fields{i});
        end
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
    
    % Set up display functions for monitoring
    fprintf('  Setting up MT-SBD...\n');
    if show_progress
        figure('Name', 'MT-SBD Progress');
        dispfun = cell(1, params.num_kernels);
        for n = 1:params.num_kernels
            dispfun{n} = @(Y, A, X, kernel_sizes, kplus) ...
                showims(data.Y_ref, data.A1{n}, data.X0_ref(:,:,n), A, X, kernel_sizes, kplus, 1);
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
        data.Y_ref, kernel_sizes_ref, sbd_params, dispfun, data.A1, initial_iteration, maxIT);
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
    
    % Store results in data struct
    data.A_ref = A_ref;
    data.X_ref = X_ref;
    data.b_ref = bout;
    
    % Store results struct (comprehensive results for analysis)
    results = struct();
    results.A = A_ref;
    results.X = X_ref;
    results.b = bout;
    results.extras = extras;
    results.params = sbd_params;
    results.mtsbd_time = mtsbd_time;
    results.final_metrics = final_metrics;
    results.final_kernel_quality = final_kernel_quality;
    
    % Update params struct with MT-SBD settings
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
    
    fprintf('  Reference slice MT-SBD complete.\n');
end


