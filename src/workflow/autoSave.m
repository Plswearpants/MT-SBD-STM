function [meta] = autoSave(log, data, params, meta, varargin)
%AUTOSAVE Intelligent automatic saving based on workflow phase
%
%   Automatically detects the current workflow phase and saves appropriately:
%   - Pre-run phase: saves dataset (auto.mat or manualXX.mat)
%   - Post-run phase: saves run results (run_<dataset>_XX/)
%
%   [meta] = autoSave(log, data, params, meta, ...)
%
%   INPUTS:
%       log                 - Log struct with .path and .file fields
%       data                - Data struct (content determines phase)
%       params              - Parameter struct (should be hierarchical for storage)
%       meta                - Metadata struct from createProjectStructure
%
%       OPTIONAL (Name-Value pairs):
%       'phase'             - Force specific phase: 'dataset' or 'run' 
%                             (default: auto-detect from data struct)
%       'compression'       - Use compression (-v7.3 for large files) (default: false)
%       'verbose'           - Show detailed output (default: true)
%
%   OUTPUTS:
%       meta                - Updated metadata struct with save information
%
%   PHASE DETECTION LOGIC:
%       Pre-run (dataset):  data.synGen exists, no mcsbd_* fields
%       Post-run (run):     data.slice or data.block exists (or legacy mcsbd_slice/mcsbd_block)
%
%   EXAMPLES:
%       % Auto-detect and save (after data generation)
%       meta = autoSave(log, data, params, meta);
%       % → Saves as dataset (auto.mat or manualXX.mat)
%
%       % Auto-detect and save (after algorithm run)
%       meta = autoSave(log, data, params, meta);
%       % → Saves as run (run_<dataset>_XX/)
%
%       % Force specific phase
%       meta = autoSave(log, data, params, meta, 'phase', 'dataset');
%
%   See also: saveDataset, saveRun, createProjectStructure

    % Parse input arguments
    p = inputParser;
    addRequired(p, 'log', @isstruct);
    addRequired(p, 'data', @isstruct);
    addRequired(p, 'params', @isstruct);
    addRequired(p, 'meta', @isstruct);
    addParameter(p, 'phase', '', @(x) isempty(x) || ismember(x, {'dataset', 'run'}));
    addParameter(p, 'compression', false, @islogical);
    addParameter(p, 'verbose', true, @islogical);
    parse(p, log, data, params, meta, varargin{:});
    
    % Extract parameters
    phase_forced = p.Results.phase;
    compression = p.Results.compression;
    verbose = p.Results.verbose;
    
    % Validate inputs
    if ~isfield(log, 'path') || ~isfield(log, 'file')
        error('log struct must contain .path and .file fields');
    end
    
    % Detect phase if not forced
    if isempty(phase_forced)
        phase = detectPhase(data);
        if verbose
            fprintf('  Auto-detected phase: %s\n', phase);
        end
    else
        phase = phase_forced;
        if verbose
            fprintf('  Using forced phase: %s\n', phase);
        end
    end
    
    % Call appropriate save function based on phase
    switch phase
        case 'dataset'
            % Pre-run phase: save dataset
            if verbose
                fprintf('  Saving dataset (pre-run phase)...\n');
            end
            meta = saveDataset(log, data, params, meta, ...
                'compression', compression);
            
        case 'run'
            % Post-run phase: save run
            if verbose
                fprintf('  Saving run (post-run phase)...\n');
            end
            meta = saveRun(log, data, params, meta, ...
                'compression', compression);
            
        otherwise
            error('Unknown phase: %s', phase);
    end
    
    if verbose
        fprintf('  Auto-save complete.\n');
    end
end

%% Helper function: Detect workflow phase from data struct
function phase = detectPhase(data)
    %DETECTPHASE Determine workflow phase from data struct contents
    %
    %   Phase detection logic:
    %   1. Check for post-run indicators (mcsbd_slice or mcsbd_block)
    %   2. Check for pre-run indicators (synGen with initialization)
    %   3. Error if ambiguous or unrecognized
    
    % Check for post-run phase
    has_mcsbd_slice = isfield(data, 'mcsbd_slice');
    has_mcsbd_block = isfield(data, 'mcsbd_block');
    
    % Check for pre-run phase
    has_synGen = isfield(data, 'synGen');
    has_initialization = (isfield(data, 'slice') && isfield(data.slice, 'A_init')) || (isfield(data, 'initialization') && isfield(data.initialization, 'A_init'));
    
    % Determine phase
    if has_mcsbd_slice || has_mcsbd_block
        % Post-run phase: algorithm has been executed
        phase = 'run';
        
    elseif has_initialization
        % Pre-run phase: data generated + kernels initialized
        phase = 'dataset';
        
    elseif has_synGen
        % Ambiguous: data generated but no kernel initialization yet
        error(['Cannot auto-detect phase: data.synGen exists but no initialization found. ', ...
               'Run kernel initialization first or force phase with ''phase'' parameter.']);
        
    else
        % Unrecognized data structure
        error(['Cannot auto-detect phase: unrecognized data struct. ', ...
               'Expected data.synGen (pre-run) or data.mcsbd_* (post-run).']);
    end
end

