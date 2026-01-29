function params_out = organizeParams(params_in, mode)
%ORGANIZEPARAMS Organize parameters between flat and hierarchical structures
%
%   Handles conversion between flat parameter structure (used internally by
%   functions) and hierarchical parameter structure (used for storage and
%   organization).
%
%   params_out = organizeParams(params_in, mode)
%
%   INPUTS:
%       params_in       - Input params struct (flat or hierarchical)
%       mode            - 'write' or 'extract'
%                         'write': flat → hierarchical (for saving)
%                         'extract': hierarchical → flat (for loading/using)
%
%   OUTPUTS:
%       params_out      - Organized params struct
%
%   HIERARCHICAL STRUCTURE:
%       params
%       ├─ synGen              % Synthetic data generation parameters
%       │  ├─ SNR, N_obs, observation_resolution, defect_density
%       │  ├─ num_slices, num_kernels, ref_slice
%       │  ├─ kernel_sizes, LDoS_path, normalization_type
%       │  └─ vis_generation, selected_indices, cutoff_M
%       │
%       ├─ initialization      % Kernel initialization parameters
%       │  ├─ init_method      % 'auto' or 'manualXX'
%       │  ├─ init_window      % Window function {type, sigma}
%       │  ├─ kernel_sizes     % Initialization kernel sizes
%       │  ├─ kernel_selection_type % 'selected' or 'random' (manual only)
%       │  └─ kernel_sizes_ref % Reference slice kernel sizes (manual only)
%       │  Note: kernel_centers stored in data.initialization, not params
%       │
%       ├─ mcsbd_slice         % Single-slice algorithm parameters
%       │  ├─ initial_iteration, maxIT, lambda1
%       │  ├─ phase2_enable, lambda2, nrefine, kplus_factor
%       │  ├─ signflip_threshold, xpos, getbias
%       │  └─ Xsolve_method
%       │
%       └─ mcsbd_block         % Block/all-slice algorithm parameters
%          ├─ maxIT_allslice, lambda_allslice
%          ├─ use_reference_init
%          └─ show_allslice_progress
%
%   FLAT STRUCTURE (internal use):
%       params.SNR, params.N_obs, params.num_kernels, ...
%       (All fields at top level for easy access in functions)
%
%   USAGE:
%       % After generation (flat → hierarchical for saving)
%       params_hierarchical = organizeParams(params_flat, 'write');
%
%       % After loading (hierarchical → flat for using)
%       params_flat = organizeParams(params_hierarchical, 'extract');
%
%   See also: generateSyntheticData, decomposeReferenceSlice

    % Validate mode
    if ~ismember(mode, {'write', 'extract'})
        error('Mode must be ''write'' or ''extract''');
    end
    
    if strcmp(mode, 'write')
        % WRITE MODE: Flat → Hierarchical
        params_out = flatToHierarchical(params_in);
    else
        % EXTRACT MODE: Hierarchical → Flat
        params_out = hierarchicalToFlat(params_in);
    end
end

%% Helper function: Flat → Hierarchical
function params_hier = flatToHierarchical(params_flat)
    params_hier = struct();
    
    % Define field mappings: namespace → field list
    synGen_fields = {'SNR', 'N_obs', 'observation_resolution', 'defect_density', ...
                     'num_slices', 'num_kernels', 'ref_slice', 'kernel_sizes', ...
                     'LDoS_path', 'normalization_type', 'vis_generation', ...
                     'selected_indices', 'cutoff_M'};
    
    initialization_fields = {'init_method', 'init_window', ...
                            'kernel_sizes', 'kernel_selection_type', 'kernel_sizes_ref'};
    
    mcsbd_slice_fields = {'initial_iteration', 'maxIT', 'lambda1', ...
                          'phase2_enable', 'lambda2', 'nrefine', 'kplus_factor', ...
                          'signflip_threshold', 'xpos', 'getbias', 'Xsolve_method'};
    
    mcsbd_block_fields = {'maxIT_allslice', 'lambda_allslice', ...
                          'use_reference_init', 'show_allslice_progress'};
    
    % Copy synGen fields
    if any(isfield(params_flat, synGen_fields))
        params_hier.synGen = struct();
        for i = 1:length(synGen_fields)
            field = synGen_fields{i};
            if isfield(params_flat, field)
                params_hier.synGen.(field) = params_flat.(field);
            end
        end
    end
    
    % Copy initialization fields (if they exist)
    if any(isfield(params_flat, initialization_fields))
        params_hier.initialization = struct();
        for i = 1:length(initialization_fields)
            field = initialization_fields{i};
            if isfield(params_flat, field)
                params_hier.initialization.(field) = params_flat.(field);
            end
        end
    end
    
    % Copy mcsbd_slice fields (if they exist)
    if any(isfield(params_flat, mcsbd_slice_fields))
        params_hier.mcsbd_slice = struct();
        for i = 1:length(mcsbd_slice_fields)
            field = mcsbd_slice_fields{i};
            if isfield(params_flat, field)
                params_hier.mcsbd_slice.(field) = params_flat.(field);
            end
        end
    end
    
    % Copy mcsbd_block fields (if they exist)
    if any(isfield(params_flat, mcsbd_block_fields))
        params_hier.mcsbd_block = struct();
        for i = 1:length(mcsbd_block_fields)
            field = mcsbd_block_fields{i};
            if isfield(params_flat, field)
                params_hier.mcsbd_block.(field) = params_flat.(field);
            end
        end
    end
    
    % Copy any other fields that don't fit into categories
    all_categorized_fields = [synGen_fields, initialization_fields, mcsbd_slice_fields, mcsbd_block_fields];
    flat_fields = fieldnames(params_flat);
    for i = 1:length(flat_fields)
        field = flat_fields{i};
        if ~ismember(field, all_categorized_fields)
            % Keep uncategorized fields at top level
            params_hier.(field) = params_flat.(field);
        end
    end
end

%% Helper function: Hierarchical → Flat
function params_flat = hierarchicalToFlat(params_hier)
    params_flat = struct();
    
    % Extract from synGen
    if isfield(params_hier, 'synGen')
        synGen_fields = fieldnames(params_hier.synGen);
        for i = 1:length(synGen_fields)
            field = synGen_fields{i};
            params_flat.(field) = params_hier.synGen.(field);
        end
    end
    
    % Extract from initialization
    if isfield(params_hier, 'initialization')
        initialization_fields = fieldnames(params_hier.initialization);
        for i = 1:length(initialization_fields)
            field = initialization_fields{i};
            params_flat.(field) = params_hier.initialization.(field);
        end
    end
    
    % Extract from mcsbd_slice
    if isfield(params_hier, 'mcsbd_slice')
        mcsbd_slice_fields = fieldnames(params_hier.mcsbd_slice);
        for i = 1:length(mcsbd_slice_fields)
            field = mcsbd_slice_fields{i};
            params_flat.(field) = params_hier.mcsbd_slice.(field);
        end
    end
    
    % Extract from mcsbd_block
    if isfield(params_hier, 'mcsbd_block')
        mcsbd_block_fields = fieldnames(params_hier.mcsbd_block);
        for i = 1:length(mcsbd_block_fields)
            field = mcsbd_block_fields{i};
            params_flat.(field) = params_hier.mcsbd_block.(field);
        end
    end
    
    % Copy any other top-level fields
    hier_fields = fieldnames(params_hier);
    for i = 1:length(hier_fields)
        field = hier_fields{i};
        if ~ismember(field, {'synGen', 'initialization', 'mcsbd_slice', 'mcsbd_block'})
            params_flat.(field) = params_hier.(field);
        end
    end
end

