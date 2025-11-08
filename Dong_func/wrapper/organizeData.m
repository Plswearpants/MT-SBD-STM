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
%       ├─ initialization      % Kernel initialization outputs
%       │  ├─ A_init           % Initialized kernels (from auto/manual init)
%       │  └─ kernel_centers   % Kernel centers from initialization (auto/manual)
%       │
%       ├─ mcsbd_slice         % Single-slice decomposition results
%       │  ├─ A, X, b          % Deconvolved kernels, activations, bias
%       │  ├─ extras           % Optimization details
%       │  ├─ mtsbd_time       % Execution time
%       │  ├─ final_metrics    % Activation quality metrics
%       │  ├─ final_kernel_quality % Kernel quality factors
%       │  ├─ kernel_centers   % Most isolated points (from IS01A, for proliferation)
%       │  ├─ A0_used          % Ground truth kernels used
%       │  ├─ defect_positions % All defect locations
%       │  ├─ most_isolated_points % Selected centers
%       │  ├─ used_most_isolated_points % Aligned centers
%       │  ├─ isolation_scores % Distance scores
%       │  ├─ num_defects      % Count per kernel
%       │  └─ offset           % Alignment offset
%       │
%       └─ mcsbd_block         % Block/all-slice algorithm results (future)
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
    
    initialization_fields = {'A_init', 'init_kernel_centers'};  
    
    mcsbd_slice_fields = {'A', 'X', 'b', 'extras', 'mtsbd_time', ...
                          'final_metrics', 'final_kernel_quality', ...
                          'kernel_centers', 'A0_used', 'defect_positions', ...
                          'most_isolated_points', 'used_most_isolated_points', ...
                          'isolation_scores', 'num_defects', 'offset'};
    
    mcsbd_block_fields = {};  % Future use - all-slice decomposition results
    
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
    
    % Copy initialization fields (if they exist)
    if any(isfield(data_flat, initialization_fields))
        data_hier.initialization = struct();
        for i = 1:length(initialization_fields)
            field = initialization_fields{i};
            if isfield(data_flat, field)
                % Special mapping: init_kernel_centers → kernel_centers in hierarchical
                if strcmp(field, 'init_kernel_centers')
                    data_hier.initialization.kernel_centers = data_flat.(field);
                else
                    data_hier.initialization.(field) = data_flat.(field);
                end
            end
        end
    end
    
    % Copy mcsbd_slice fields (if they exist)
    if any(isfield(data_flat, mcsbd_slice_fields))
        data_hier.mcsbd_slice = struct();
        for i = 1:length(mcsbd_slice_fields)
            field = mcsbd_slice_fields{i};
            if isfield(data_flat, field)
                data_hier.mcsbd_slice.(field) = data_flat.(field);
            end
        end
    end
    
    % Copy mcsbd_block fields (if they exist)
    if ~isempty(mcsbd_block_fields) && any(isfield(data_flat, mcsbd_block_fields))
        data_hier.mcsbd_block = struct();
        for i = 1:length(mcsbd_block_fields)
            field = mcsbd_block_fields{i};
            if isfield(data_flat, field)
                data_hier.mcsbd_block.(field) = data_flat.(field);
            end
        end
    end
    
    % Copy any other fields that don't fit into categories
    all_categorized_fields = [synGen_fields, initialization_fields, mcsbd_slice_fields, mcsbd_block_fields];
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
    
    % Extract from initialization
    if isfield(data_hier, 'initialization')
        initialization_fields = fieldnames(data_hier.initialization);
        for i = 1:length(initialization_fields)
            field = initialization_fields{i};
            % Special mapping: kernel_centers → init_kernel_centers in flat (to avoid conflict with mcsbd_slice.kernel_centers)
            if strcmp(field, 'kernel_centers')
                data_flat.init_kernel_centers = data_hier.initialization.(field);
            else
                data_flat.(field) = data_hier.initialization.(field);
            end
        end
    end
    
    % Extract from mcsbd_slice
    if isfield(data_hier, 'mcsbd_slice')
        mcsbd_slice_fields = fieldnames(data_hier.mcsbd_slice);
        for i = 1:length(mcsbd_slice_fields)
            field = mcsbd_slice_fields{i};
            data_flat.(field) = data_hier.mcsbd_slice.(field);
        end
    end
    
    % Extract from mcsbd_block
    if isfield(data_hier, 'mcsbd_block')
        mcsbd_block_fields = fieldnames(data_hier.mcsbd_block);
        for i = 1:length(mcsbd_block_fields)
            field = mcsbd_block_fields{i};
            data_flat.(field) = data_hier.mcsbd_block.(field);
        end
    end
    
    % Copy any other top-level fields
    hier_fields = fieldnames(data_hier);
    for i = 1:length(hier_fields)
        field = hier_fields{i};
        if ~ismember(field, {'synGen', 'initialization', 'mcsbd_slice', 'mcsbd_block'})
            data_flat.(field) = data_hier.(field);
        end
    end
end

