function data_out = organizeData(data_in, mode)
%ORGANIZEDATA Organize data between flat and hierarchical structures
%
%   Handles conversion between flat data structure (used internally by
%   functions) and hierarchical data structure (used for storage and
%   organization).
%
%   data_out = organizeData(data_in, mode)
%
%   INPUTS:
%       data_in       - Input data struct (flat or hierarchical)
%       mode          - 'write' or 'extract'
%                       'write': flat → hierarchical (for saving)
%                       'extract': hierarchical → flat (for loading/using)
%
%   OUTPUTS:
%       data_out      - Organized data struct
%
%   HIERARCHICAL STRUCTURE:
%       data
%       ├─ synGen              % Synthetic data generation outputs
%       │  ├─ Y, A0, A0_noiseless, X0
%       │  └─ Y_ref, X0_ref, A0_ref
%       │
%       ├─ slice               % Single-slice: init + decomposition (flattened)
%       │  ├─ A_init, init_kernel_centers  % From kernel init (auto/manual)
%       │  ├─ A, X, b          % Deconvolved kernels, activations, bias
%       │  ├─ extras, mtsbd_time, final_metrics, final_kernel_quality
%       │  ├─ kernel_centers   % Most isolated points (IS01A, for block init proliferation)
%       │  ├─ A0_used, defect_positions, most_isolated_points
%       │  ├─ used_most_isolated_points, isolation_scores, num_defects, offset
%       │
%       └─ block               % Block init (IB01A) + all-slice results (DA01A)
%          Block init (unified: proliferation or block_manual) and DA01A write here.
%          Block-init fields (from IB01A):
%          ├─ A_all_init           % {S×K} initialized kernels per slice
%          ├─ A_all_init_matrix    % {K×1} 3D [H×W×S] unified (default unify: max_per_kernel)
%          ├─ init_kernel_centers  % [K×2] centers (proliferation) or empty (block_manual)
%          DA01A result fields:
%          ├─ Aout_all, Xout_all, bout_all, extras_all
%
%   FLAT STRUCTURE (internal use):
%       data.Y, data.A0, data.X, data.A, ...
%       (All fields at top level for easy access in functions)
%
%   USAGE:
%       % After generation (flat → hierarchical for saving)
%       data_hierarchical = organizeData(data_flat, 'write');
%
%       % After loading (hierarchical → flat for using)
%       data_flat = organizeData(data_hierarchical, 'extract');
%
%   See also: organizeParams, generateSyntheticData, decomposeReferenceSlice

    % Validate mode
    if ~ismember(mode, {'write', 'extract'})
        error('Mode must be ''write'' or ''extract''');
    end
    
    if strcmp(mode, 'write')
        % WRITE MODE: Flat → Hierarchical
        data_out = flatToHierarchical(data_in);
    else
        % EXTRACT MODE: Hierarchical → Flat
        data_out = hierarchicalToFlat(data_in);
    end
end

%% Helper function: Flat → Hierarchical
function data_hier = flatToHierarchical(data_flat)
    data_hier = struct();
    
    % Define field mappings: namespace → field list
    synGen_fields = {'Y', 'A0', 'A0_noiseless', 'X0', ...
                     'Y_ref', 'X0_ref', 'A0_ref'};
    
    % slice = init (ref-slice) + single-slice decomposition, flattened
    slice_fields = {'A_init', 'init_kernel_centers', ...
                    'A', 'X', 'b', 'extras', 'mtsbd_time', ...
                    'final_metrics', 'final_kernel_quality', ...
                    'kernel_centers', 'A0_used', 'defect_positions', ...
                    'most_isolated_points', 'used_most_isolated_points', ...
                    'isolation_scores', 'num_defects', 'offset'};
    
    % block = block-init (IB01A) + all-slice results (DA01A). Flat names: A1_all, A1_all_matrix,
    % proliferation_kernel_centers (block init) → block.A_all_init, A_all_init_matrix, init_kernel_centers
    block_init_flat = {'A1_all', 'A1_all_matrix', 'proliferation_kernel_centers'};
    block_result_flat = {'Aout_all', 'Xout_all', 'bout_all', 'extras_all'};
    block_fields = [block_init_flat, block_result_flat];
    
    % Copy synGen fields
    if any(isfield(data_flat, synGen_fields))
        data_hier.synGen = struct();
        for i = 1:length(synGen_fields)
            field = synGen_fields{i};
            if isfield(data_flat, field)
                data_hier.synGen.(field) = data_flat.(field);
            end
        end
    end
    
    % Copy slice fields (init + single-slice decomposition)
    if any(isfield(data_flat, slice_fields))
        data_hier.slice = struct();
        for i = 1:length(slice_fields)
            field = slice_fields{i};
            if isfield(data_flat, field)
                data_hier.slice.(field) = data_flat.(field);
            end
        end
    end
    
    % Copy block fields (init + results); block-init flat names map to block hier names
    if any(isfield(data_flat, block_fields))
        data_hier.block = struct();
        for i = 1:length(block_init_flat)
            field = block_init_flat{i};
            if isfield(data_flat, field)
                if strcmp(field, 'proliferation_kernel_centers')
                    data_hier.block.init_kernel_centers = data_flat.(field);
                elseif strcmp(field, 'A1_all')
                    data_hier.block.A_all_init = data_flat.(field);
                elseif strcmp(field, 'A1_all_matrix')
                    data_hier.block.A_all_init_matrix = data_flat.(field);
                else
                    data_hier.block.(field) = data_flat.(field);
                end
            end
        end
        for i = 1:length(block_result_flat)
            field = block_result_flat{i};
            if isfield(data_flat, field)
                data_hier.block.(field) = data_flat.(field);
            end
        end
    end
    
    % Copy any other fields that don't fit into categories
    all_categorized_fields = [synGen_fields, slice_fields, block_fields];
    flat_fields = fieldnames(data_flat);
    for i = 1:length(flat_fields)
        field = flat_fields{i};
        if ~ismember(field, all_categorized_fields)
            % Keep uncategorized fields at top level
            data_hier.(field) = data_flat.(field);
        end
    end
end

%% Helper function: Hierarchical → Flat
function data_flat = hierarchicalToFlat(data_hier)
    data_flat = struct();
    
    % Extract from synGen
    if isfield(data_hier, 'synGen')
        synGen_fields = fieldnames(data_hier.synGen);
        for i = 1:length(synGen_fields)
            field = synGen_fields{i};
            data_flat.(field) = data_hier.synGen.(field);
        end
    end
    
    % Extract from slice (init + single-slice decomposition)
    if isfield(data_hier, 'slice')
        slice_fields = fieldnames(data_hier.slice);
        for i = 1:length(slice_fields)
            field = slice_fields{i};
            data_flat.(field) = data_hier.slice.(field);
        end
    end
    
    % Backward compatibility: old files may have initialization / mcsbd_slice
    if isfield(data_hier, 'initialization')
        init_f = fieldnames(data_hier.initialization);
        for i = 1:length(init_f)
            f = init_f{i};
            if strcmp(f, 'kernel_centers')
                data_flat.init_kernel_centers = data_hier.initialization.(f);
            else
                data_flat.(f) = data_hier.initialization.(f);
            end
        end
    end
    if isfield(data_hier, 'mcsbd_slice')
        slice_f = fieldnames(data_hier.mcsbd_slice);
        for i = 1:length(slice_f)
            data_flat.(slice_f{i}) = data_hier.mcsbd_slice.(slice_f{i});
        end
    end
    
    % Extract from proliferation
    if isfield(data_hier, 'proliferation')
        prolif_f = fieldnames(data_hier.proliferation);
        for i = 1:length(prolif_f)
            field = prolif_f{i};
            if strcmp(field, 'kernel_centers')
                data_flat.proliferation_kernel_centers = data_hier.proliferation.(field);
            else
                data_flat.(field) = data_hier.proliferation.(field);
            end
        end
    end
    
    % Extract from block (block-init names → flat names for downstream)
    if isfield(data_hier, 'block')
        block_f = fieldnames(data_hier.block);
        for i = 1:length(block_f)
            field = block_f{i};
            if strcmp(field, 'A_all_init')
                data_flat.A1_all = data_hier.block.(field);
            elseif strcmp(field, 'A_all_init_matrix')
                data_flat.A1_all_matrix = data_hier.block.(field);
            elseif strcmp(field, 'init_kernel_centers')
                data_flat.proliferation_kernel_centers = data_hier.block.(field);
            else
                data_flat.(field) = data_hier.block.(field);
            end
        end
    end
    if isfield(data_hier, 'mcsbd_block')
        block_f = fieldnames(data_hier.mcsbd_block);
        for i = 1:length(block_f)
            data_flat.(block_f{i}) = data_hier.mcsbd_block.(block_f{i});
        end
    end
    
    % Copy any other top-level fields
    hier_fields = fieldnames(data_hier);
    for i = 1:length(hier_fields)
        field = hier_fields{i};
        if ~ismember(field, {'synGen', 'slice', 'block', 'initialization', 'mcsbd_slice', 'mcsbd_block', 'proliferation'})
            data_flat.(field) = data_hier.(field);
        end
    end
end

