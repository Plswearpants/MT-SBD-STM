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
%       ├─ slice               % Single-slice: init + decomposition (flattened)
%       │  ├─ init_method, init_window, kernel_sizes_ref, kernel_selection_type
%       │  ├─ initial_iteration, maxIT, lambda1
%       │  ├─ phase2_enable, lambda2, nrefine, kplus_factor
%       │  ├─ signflip_threshold, xpos, getbias, Xsolve_method
%       │  └─ use_xinit, show_progress
%       │  Note: init_kernel_centers stored in data.slice, not params
%       │
%       └─ block               % Block init (IB01A) + all-slice algo (DA01A)
%          Block-init params (IB01A): method, A1_matrix_unify_size, window_type, ...
%          DA01A params: maxIT_allslice, use_reference_init, show_allslice_progress
%          ├─ block_init_method  % 'proliferation' | 'block_manual' (unified IB01A)
%          ├─ A1_matrix_unify_size % 'max_per_kernel' (default) | 'ref_slice'
%          ├─ window_type_proliferation, interactive_size_adjust, use_matrix_format
%          ├─ maxIT_allslice, lambda_allslice
%          ├─ use_reference_init, show_allslice_progress
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
    
    % slice = initialization (ref-slice init) + mcsbd_slice (single-slice algo), flattened
    slice_fields = {'init_method', 'init_window', 'kernel_sizes', 'kernel_selection_type', 'kernel_sizes_ref', ...
                   'initial_iteration', 'maxIT', 'lambda1', 'phase2_enable', 'lambda2', 'nrefine', 'kplus_factor', ...
                   'signflip_threshold', 'xpos', 'getbias', 'Xsolve_method', 'use_xinit', 'show_progress'};
    
    block_fields = {'block_init_method', 'A1_matrix_unify_size', ...
                    'window_type_proliferation', 'interactive_size_adjust', 'use_matrix_format', ...
                    'initial_iteration', 'maxIT_allslice', 'lambda1', 'phase2_enable', 'lambda2', ...
                    'nrefine', 'kplus_factor', 'signflip_threshold', 'xpos', 'getbias', ...
                    'Xsolve_method', 'use_xinit', ...
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
    
    % Copy slice fields (init + single-slice algo), flattened
    if any(isfield(params_flat, slice_fields))
        params_hier.slice = struct();
        for i = 1:length(slice_fields)
            field = slice_fields{i};
            if isfield(params_flat, field)
                params_hier.slice.(field) = params_flat.(field);
            end
        end
    end
    
    % Copy block fields (if they exist)
    if any(isfield(params_flat, block_fields))
        params_hier.block = struct();
        for i = 1:length(block_fields)
            field = block_fields{i};
            if isfield(params_flat, field)
                params_hier.block.(field) = params_flat.(field);
            end
        end
    end
    
    % Copy any other fields that don't fit into categories
    all_categorized_fields = [synGen_fields, slice_fields, block_fields];
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
    
    % Extract from slice (init + single-slice algo)
    if isfield(params_hier, 'slice')
        slice_fields = fieldnames(params_hier.slice);
        for i = 1:length(slice_fields)
            field = slice_fields{i};
            params_flat.(field) = params_hier.slice.(field);
        end
    end
    
    % Backward compatibility: old files may have initialization / mcsbd_slice
    if isfield(params_hier, 'initialization')
        init_f = fieldnames(params_hier.initialization);
        for i = 1:length(init_f)
            params_flat.(init_f{i}) = params_hier.initialization.(init_f{i});
        end
    end
    if isfield(params_hier, 'mcsbd_slice')
        slice_f = fieldnames(params_hier.mcsbd_slice);
        for i = 1:length(slice_f)
            params_flat.(slice_f{i}) = params_hier.mcsbd_slice.(slice_f{i});
        end
    end
    
    % Extract from block
    if isfield(params_hier, 'block')
        block_fields = fieldnames(params_hier.block);
        for i = 1:length(block_fields)
            field = block_fields{i};
            params_flat.(field) = params_hier.block.(field);
        end
    end
    if isfield(params_hier, 'mcsbd_block')
        block_f = fieldnames(params_hier.mcsbd_block);
        for i = 1:length(block_f)
            params_flat.(block_f{i}) = params_hier.mcsbd_block.(block_f{i});
        end
    end
    
    % Copy any other top-level fields
    hier_fields = fieldnames(params_hier);
    for i = 1:length(hier_fields)
        field = hier_fields{i};
        if ~ismember(field, {'synGen', 'slice', 'block', 'initialization', 'mcsbd_slice', 'mcsbd_block'})
            params_flat.(field) = params_hier.(field);
        end
    end
end

