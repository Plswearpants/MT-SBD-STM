function [data, params] = generateSyntheticData(varargin)
%GENERATESYNTHETICDATA Wrapper function for synthetic data generation
%
%   Complete workflow for generating synthetic STM data including:
%   - 1. Data generation from LDoS simulations
%   - 2. Normalization and projection
%   - 3. Reference slice selection
%   - 4. Data extraction and organization
%
%   INPUTS (Name-Value pairs):
%       REQUIRED (will prompt if not provided):
%       'SNR'                   - Signal-to-noise ratio
%       'N_obs'                 - Observation lattice size
%       'observation_resolution' - Pixels per lattice site
%       'defect_density'        - Surface defect density (0-1)
%       'num_slices'            - Number of energy slices
%       'LDoS_path'             - Path to LDoS data file
%
%       OPTIONAL (have defaults):
%       'vis_generation'        - Show intermediate steps (default: false)
%       'normalization_type'    - 'dynamic' or 'static' (default: 'dynamic')
%       'ref_slice'             - Reference slice number (default: [] for interactive)
%
%   OUTPUTS:
%       data            - Structure containing all data arrays:
%           .Y              - Normalized 3D observation [H x W x num_slices]
%           .A0             - Noisy kernels cell array {num_slices x num_kernels}
%           .A0_noiseless   - Noiseless kernels cell array {num_slices x num_kernels}
%           .X0             - Ground truth activations [H x W x num_kernels]
%           .Y_ref          - Reference slice observation [H x W]
%           .X0_ref         - Reference slice activations (same as X0)
%           .A0_ref         - Reference slice kernels {1 x num_kernels}
%
%       params          - Structure containing all parameters (1D scalars/strings):
%           .SNR, .N_obs, .observation_resolution, .defect_density
%           .num_slices, .num_kernels, .ref_slice, etc.
%
%   EXAMPLE:
%       [data, params] = generateSyntheticData(...
%           'SNR', 5, 'N_obs', 50, 'num_slices', 2, ...
%           'LDoS_path', 'LDoS_single_defects_self=0.6_save.mat');
%       Y = data.Y;
%       A0 = data.A0;
%
%   See also: properGen_full, normalizeBackgroundToZeroMean3D, proj2oblique

    % Parse input arguments
    p = inputParser;
    % Required parameters (user must provide or will be prompted)
    addParameter(p, 'SNR', [], @(x) isempty(x) || isnumeric(x));
    addParameter(p, 'N_obs', [], @(x) isempty(x) || isnumeric(x));
    addParameter(p, 'observation_resolution', [], @(x) isempty(x) || isnumeric(x));
    addParameter(p, 'defect_density', [], @(x) isempty(x) || isnumeric(x));
    addParameter(p, 'num_slices', [], @(x) isempty(x) || isnumeric(x));
    addParameter(p, 'LDoS_path', '', @ischar);
    % Optional parameters (have defaults)
    addParameter(p, 'vis_generation', false, @islogical);
    addParameter(p, 'normalization_type', 'dynamic', @ischar);
    addParameter(p, 'ref_slice', [], @isnumeric);
    parse(p, varargin{:});
    
    % Store all input parameters in params struct
    params = struct();
    
    % Handle required parameters - prompt if not provided
    if isempty(p.Results.SNR)
        params.SNR = input('Enter SNR (Signal-to-Noise Ratio): ');
        if isempty(params.SNR) || ~isnumeric(params.SNR) || params.SNR <= 0
            error('SNR must be a positive number');
        end
    else
        params.SNR = p.Results.SNR;
    end
    
    if isempty(p.Results.N_obs)
        params.N_obs = input('Enter N_obs (Observation lattice size): ');
        if isempty(params.N_obs) || ~isnumeric(params.N_obs) || params.N_obs <= 0
            error('N_obs must be a positive integer');
        end
    else
        params.N_obs = p.Results.N_obs;
    end
    
    if isempty(p.Results.observation_resolution)
        params.observation_resolution = input('Enter observation_resolution (Pixels per lattice site): ');
        if isempty(params.observation_resolution) || ~isnumeric(params.observation_resolution) || params.observation_resolution <= 0
            error('observation_resolution must be a positive integer');
        end
    else
        params.observation_resolution = p.Results.observation_resolution;
    end
    
    if isempty(p.Results.defect_density)
        params.defect_density = input('Enter defect_density (0-1): ');
        if isempty(params.defect_density) || ~isnumeric(params.defect_density) || params.defect_density <= 0 || params.defect_density >= 1
            error('defect_density must be between 0 and 1');
        end
    else
        params.defect_density = p.Results.defect_density;
    end
    
    if isempty(p.Results.num_slices)
        params.num_slices = input('Enter num_slices (Number of energy slices): ');
        if isempty(params.num_slices) || ~isnumeric(params.num_slices) || params.num_slices <= 0
            error('num_slices must be a positive integer');
        end
    else
        params.num_slices = p.Results.num_slices;
    end
    
    if isempty(p.Results.LDoS_path)
        params.LDoS_path = input('Enter LDoS_path (Path to LDoS data file): ', 's');
        if isempty(params.LDoS_path)
            error('LDoS_path is required');
        end
    else
        params.LDoS_path = p.Results.LDoS_path;
    end
    
    % Optional parameters with defaults
    params.vis_generation = p.Results.vis_generation;
    params.normalization_type = p.Results.normalization_type;
    params.ref_slice_input = p.Results.ref_slice;  % Store user input
    
    % Generate synthetic data
    fprintf('  Generating synthetic data...\n');
    [Y, A0_noiseless, X0, gen_params] = properGen_full(...
        params.SNR, params.N_obs, params.observation_resolution, params.defect_density, ...
        'LDoS_path', params.LDoS_path, ...
        'num_slices', params.num_slices, ...
        'vis', params.vis_generation);
    
    % Merge generation parameters into params (scalars/1D only)
    params.num_kernels = gen_params.num_kernels;
    params.num_slices = gen_params.num_slices;
    if isfield(gen_params, 'selected_indices')
        params.selected_indices = gen_params.selected_indices;
    end
    if isfield(gen_params, 'cutoff_M')
        params.cutoff_M = gen_params.cutoff_M;
    end
    
    % Standardize kernel sizes structure (store in params as metadata)
    params.kernel_sizes = zeros(params.num_slices, params.num_kernels, 2);
    for s = 1:params.num_slices
        for k = 1:params.num_kernels
            params.kernel_sizes(s,k,:) = size(gen_params.A0{s,k});
        end
    end
    
    % Normalize observation
    fprintf('  Normalizing observation...\n');
    Y = normalizeBackgroundToZeroMean3D(Y, params.normalization_type);
    Y = proj2oblique(Y);
    
    fprintf('  Generated %d kernels across %d energy slices.\n', params.num_kernels, params.num_slices);
    fprintf('  Observation size: %dx%dx%d pixels\n', size(Y,1), size(Y,2), size(Y,3));
    
    % Display and select reference slice
    if isempty(params.ref_slice_input)
        % Interactive selection
        figure('Name', 'Generated 3D Data');
        d3gridDisplay(Y, params.normalization_type);
        title('Select a reference slice for initial decomposition');
        
        params.ref_slice = input('Enter reference slice number: ');
        
        % Validate reference slice
        if isempty(params.ref_slice) || ~isnumeric(params.ref_slice) || ...
           params.ref_slice < 1 || params.ref_slice > size(Y, 3)
            error('Invalid reference slice. Please enter a number between 1 and %d.', size(Y, 3));
        end
    else
        % Use provided ref_slice
        params.ref_slice = params.ref_slice_input;
        % Validate provided ref_slice
        if params.ref_slice < 1 || params.ref_slice > size(Y, 3)
            error('Invalid reference slice %d. Must be between 1 and %d.', params.ref_slice, size(Y, 3));
        end
    end
    
    fprintf('  Using slice %d as reference slice.\n', params.ref_slice);
    
    % Extract reference slice data
    Y_ref = Y(:,:,params.ref_slice);
    X0_ref = X0;
    A0_ref = cell(1, params.num_kernels);
    for i = 1:params.num_kernels
        A0_ref{i} = gen_params.A0{params.ref_slice,i};
    end
    
    % Display reference slice
    figure('Name', 'Reference Slice');
    imagesc(Y_ref);
    colorbar;
    title(sprintf('Reference Slice %d', params.ref_slice));
    axis square;
    colormap(gray);
    
    % Organize outputs into data and params structs
    % DATA struct: all multi-dimensional arrays
    data = struct();
    data.Y = Y;
    data.A0 = gen_params.A0;
    data.A0_noiseless = A0_noiseless;
    data.X0 = X0;
    data.Y_ref = Y_ref;
    data.X0_ref = X0_ref;
    data.A0_ref = A0_ref;
    
    % PARAMS struct already contains all parameters (1D scalars/strings)
    % Remove temporary field
    params = rmfield(params, 'ref_slice_input');
    
    fprintf('  Data generation complete.\n');
end

