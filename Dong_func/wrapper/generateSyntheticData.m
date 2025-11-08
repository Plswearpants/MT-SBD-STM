function [data, params] = generateSyntheticData(log, params)
%GENERATESYNTHETICDATA Wrapper function for synthetic data generation
%   Now includes internal logging
%
%   Complete workflow for generating synthetic STM data including:
%   - 1. Data generation from LDoS simulations
%   - 2. Normalization and projection
%   - 3. Reference slice selection
%   - 4. Data extraction and organization
%
%   INPUTS:
%       log                 - Log struct with .path and .file fields
%       params              - Parameter struct (hierarchical) with fields:
%
%       REQUIRED (in params.synGen):
%       .SNR                - Signal-to-noise ratio
%       .N_obs              - Observation lattice size
%       .observation_resolution - Pixels per lattice site
%       .defect_density     - Surface defect density (0-1)
%       .num_slices         - Number of energy slices
%       .LDoS_path          - Path to LDoS data file
%
%       OPTIONAL (have defaults):
%       .vis_generation     - Show intermediate steps (default: false)
%       .normalization_type - 'dynamic' or 'static' (default: 'dynamic')
%       .ref_slice          - Reference slice number (default: [] for interactive)
%
%   OUTPUTS:
%       data            - Structure containing synthetic generation data:
%           .synGen.Y              - Normalized 3D observation [H x W x num_slices]
%           .synGen.A0             - Noisy kernels cell array {num_slices x num_kernels}
%           .synGen.A0_noiseless   - Noiseless kernels cell array {num_slices x num_kernels}
%           .synGen.X0             - Ground truth activations [H x W x num_kernels]
%           .synGen.Y_ref          - Reference slice observation [H x W]
%           .synGen.X0_ref         - Reference slice activations (same as X0)
%           .synGen.A0_ref         - Reference slice kernels {1 x num_kernels}
%
%       params          - Structure containing all parameters (hierarchical for storage):
%           .synGen.SNR, .N_obs, .observation_resolution, .defect_density
%           .synGen.num_slices, .num_kernels, .ref_slice, .kernel_sizes, etc.
%           Note: Organized into hierarchical structure via organizeParams()
%
%   EXAMPLE:
%       % Set parameters in hierarchical structure
%       params.synGen.SNR = 5;
%       params.synGen.N_obs = 50;
%       params.synGen.num_slices = 2;
%       params.synGen.LDoS_path = 'LDoS_single_defects_self=0.6_save.mat';
%       [data, params] = generateSyntheticData(log, params);
%       Y = data.synGen.Y;
%       A0 = data.synGen.A0;
%
%   See also: properGen_full, normalizeBackgroundToZeroMean3D, proj2oblique

    % Validate required inputs
    if ~isstruct(log) || ~isfield(log, 'path') || ~isfield(log, 'file')
        error('log must be a struct with .path and .file fields');
    end
    if ~isstruct(params)
        error('params must be a struct');
    end
    
    % Extract params from hierarchical structure (if needed)
    if isfield(params, 'synGen')
        params = organizeParams(params, 'extract');
    end
    
    % Read required parameters from params struct (all required, no defaults)
    required_fields = {'SNR', 'N_obs', 'observation_resolution', 'defect_density', ...
                       'num_slices', 'LDoS_path'};
    for i = 1:length(required_fields)
        field = required_fields{i};
        if ~isfield(params, field)
            error('Parameter "%s" is required in params struct.', field);
        end
        if isempty(params.(field))
            error('Parameter "%s" cannot be empty.', field);
        end
    end
    
    % Validate required parameters
    if ~isnumeric(params.SNR) || params.SNR <= 0
        error('SNR must be a positive number');
    end
    if ~isnumeric(params.N_obs) || params.N_obs <= 0
        error('N_obs must be a positive integer');
    end
    if ~isnumeric(params.observation_resolution) || params.observation_resolution <= 0
        error('observation_resolution must be a positive integer');
    end
    if ~isnumeric(params.defect_density) || params.defect_density <= 0 || params.defect_density >= 1
        error('defect_density must be between 0 and 1');
    end
    if ~isnumeric(params.num_slices) || params.num_slices <= 0
        error('num_slices must be a positive integer');
    end
    if ~ischar(params.LDoS_path) && ~isstring(params.LDoS_path)
        error('LDoS_path must be a string');
    end
    
    % Read optional parameters with defaults
    if isfield(params, 'vis_generation'), vis_generation = params.vis_generation; else, vis_generation = false; end
    if isfield(params, 'normalization_type'), normalization_type = params.normalization_type; else, normalization_type = 'dynamic'; end
    if isfield(params, 'ref_slice'), ref_slice_input = params.ref_slice; else, ref_slice_input = []; end
    
    % Store in params struct for internal use (flat structure)
    params.vis_generation = vis_generation;
    params.normalization_type = normalization_type;
    params.ref_slice_input = ref_slice_input;  % Store user input
    
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
    % DATA struct: Store flat internally (will convert to hierarchical at end)
    data = struct();
    data.Y = Y;
    data.A0 = gen_params.A0;
    data.A0_noiseless = A0_noiseless;
    data.X0 = X0;
    data.Y_ref = Y_ref;
    data.X0_ref = X0_ref;
    data.A0_ref = A0_ref;
    
    % PARAMS struct: Remove temporary field first
    params = rmfield(params, 'ref_slice_input');
    
    % LOG: Generation details
    LOGcomment = sprintf("SNR=%g, N_obs=%d, resolution=%d, density=%g, slices_chosen= %d, slices=%d", ...
        params.SNR, params.N_obs, params.observation_resolution, params.defect_density, params.num_slices);
    LOGcomment = logUsedBlocks(log.path, log.file, "GD01A", LOGcomment, 0);
    
    LOGcomment = sprintf("Generated Y: %dx%dx%d, Kernels: %d, Ref slice: %d", ...
        size(data.Y,1), size(data.Y,2), size(data.Y,3), params.num_kernels, params.ref_slice);
    LOGcomment = logUsedBlocks(log.path, log.file, "  ^  ", LOGcomment, 0);
    
    % Convert data and params to hierarchical structure for storage
    data = organizeData(data, 'write');
    params = organizeParams(params, 'write');
    
    fprintf('  Data generation complete.\n');
end

