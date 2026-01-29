function [log, data, params, meta] = loadSequential(varargin)
%LOADSEQUENTIAL Sequential loading with hierarchical merging of dataset and run files
%
%   Loads workspace files with automatic hierarchical merging. If loading a
%   run file, automatically loads and merges the parent dataset file.
%
%   [log, data, params, meta] = loadSequential(...)
%
%   INPUTS (all optional Name-Value pairs):
%       'workspace_file'    - Full path to .mat file to load (default: UI selection)
%       'current_log'       - Current log struct (for memory protection)
%       'current_data'      - Current data struct (for save prompt)
%       'current_params'    - Current params struct (for save prompt)
%
%   OUTPUTS:
%       log                 - Merged log struct (most recent takes precedence)
%       data                - Merged data struct with all substructs
%       params              - Merged params struct (later overwrites earlier)
%       meta                - Metadata struct with load history
%
%   LOADING LOGIC:
%       1. Load the selected file
%       2. Detect file type (dataset vs run) based on location
%       3. If run file: load parent dataset and merge data structs
%       4. Merge all data substructs into single struct
%       5. Use most recent log, merge params (later overwrites earlier)
%
%   DATA MERGING:
%       - Dataset file contains: data.synGen
%       - Run file contains: data.mcsbd_slice or data.mcsbd_block
%       - Merged result: data.synGen + data.mcsbd_slice/block
%
%   DESCRIPTION:
%       This function provides intelligent loading for the hierarchical
%       project structure. It automatically:
%       - Detects if loading a run file and loads parent dataset
%       - Merges data structs to provide complete workspace
%       - Handles memory protection for existing workspaces
%       - Tracks load history in meta struct
%
%   EXAMPLE:
%       % Load with UI selection
%       [log, data, params, meta] = loadSequential();
%
%       % Load specific run file (auto-loads parent dataset)
%       [log, data, params, meta] = loadSequential('workspace_file', ...
%           'C:/projects/synthetic_20251031_120000/run01/run01.mat');
%
%       % Load with existing workspace protection
%       [log, data, params, meta] = loadSequential('current_log', log, ...
%           'current_data', data, 'current_params', params);
%
%   See also: saveDataset, saveRun, createProjectStructure, loadWorkspace

    % Parse input arguments
    p = inputParser;
    addParameter(p, 'workspace_file', '', @ischar);
    addParameter(p, 'current_log', struct(), @isstruct);
    addParameter(p, 'current_data', struct(), @isstruct);
    addParameter(p, 'current_params', struct(), @isstruct);
    parse(p, varargin{:});
    
    % Extract parameters
    workspace_file = p.Results.workspace_file;
    current_log = p.Results.current_log;
    current_data = p.Results.current_data;
    current_params = p.Results.current_params;
    
    fprintf('========================================================================\n');
    fprintf('  LOADING WORKSPACE (Sequential Hierarchical Loading)\n');
    fprintf('========================================================================\n');
    
    % Check if current workspace exists
    has_current_workspace = isstruct(current_log) && ...
                           isfield(current_log, 'file') && ...
                           ~isempty(current_log.file);
    
    if has_current_workspace
        fprintf('Current workspace detected in memory.\n');
        fprintf('Loading will overwrite current workspace in memory.\n');
        fprintf('WARNING: Current workspace is NOT automatically saved.\n');
        fprintf('If you need to save it, cancel this operation and use saveDataset/saveRun first.\n\n');
        
        % Prompt: Delete or rename base name in memory?
        fprintf('How do you want to handle the current workspace in memory?\n');
        fprintf('  1. Delete (clear the base name)\n');
        fprintf('  2. Rename (add suffix to base name)\n');
        fprintf('  3. Cancel load operation\n');
        
        choice = input('Enter choice (1-3): ', 's');
        
        switch choice
            case '1'
                % Delete - just clear it (will be replaced by loaded workspace)
                clear current_log;
                clear current_data;
                clear current_params;
                fprintf('Current workspace deleted.\n');
                
            case '2'
                % Rename - create new variables with suffix in base workspace
                suffix = input('Enter suffix to add to variable names: ', 's');
                if isempty(suffix)
                    fprintf('No suffix provided. default suffix "_hist" will be added.\n');
                    suffix = '_hist';
                end
                
                % Create renamed variables in base workspace
                new_log_name = ['log' suffix];
                new_data_name = ['data' suffix];
                new_params_name = ['params' suffix];
                
                assignin('base', new_log_name, current_log);
                assignin('base', new_data_name, current_data);
                assignin('base', new_params_name, current_params);
                
                fprintf('Current workspace variables renamed:\n');
                fprintf('  log → %s\n', new_log_name);
                fprintf('  data → %s\n', new_data_name);
                fprintf('  params → %s\n', new_params_name);
                fprintf('Original variables (log, data, params) will be replaced by loaded workspace.\n');
                
            case '3'
                error('Load cancelled by user.');
                
            otherwise
                error('Invalid choice. Load cancelled.');
        end
        
        fprintf('Proceeding with load...\n\n');
    else
        fprintf('No current workspace detected. Proceeding with fresh load.\n\n');
    end
    
    % Select workspace file if not provided
    if isempty(workspace_file)
        fprintf('Select workspace file to load (.mat file)...\n');
        [file, path] = uigetfile('*.mat', 'Select workspace file');
        if isequal(file, 0)
            error('No file selected. Load cancelled.');
        end
        workspace_file = fullfile(path, file);
    end
    
    fprintf('Loading from: %s\n', workspace_file);
    
    % Validate file exists
    if ~exist(workspace_file, 'file')
        error('Workspace file not found: %s', workspace_file);
    end
    
    % Parse file path
    [file_dir, file_name, ~] = fileparts(workspace_file);
    [parent_dir, current_folder] = fileparts(file_dir);
    
    % Detect file type based on location and name
    is_run_file = startsWith(current_folder, 'run') && startsWith(file_name, 'run');
    is_dataset_file = ~is_run_file && (strcmp(file_name, 'auto') || startsWith(file_name, 'manual'));
    
    % Initialize load history
    load_history = struct();
    load_history.files_loaded = {};
    load_history.load_order = {};
    
    % Load the selected file
    fprintf('  Loading selected file...\n');
    loaded = load(workspace_file);
    
    if ~isfield(loaded, 'log') || ~isfield(loaded, 'data') || ~isfield(loaded, 'params')
        error('Invalid workspace file: must contain log, data, and params');
    end
    
    % Store loaded data
    log = loaded.log;
    data = loaded.data;
    params = loaded.params;
    load_history.files_loaded{end+1} = workspace_file;
    load_history.load_order{end+1} = file_name;
    
    fprintf('    Loaded: %s\n', file_name);
    
    % If this is a run file, load parent dataset
    if is_run_file
        fprintf('  Detected run file. Loading parent dataset...\n');
        
        % Find dataset file in parent directory
        dataset_files = [dir(fullfile(parent_dir, 'auto.mat')); ...
                        dir(fullfile(parent_dir, 'manual*.mat'))];
        
        if isempty(dataset_files)
            warning('No dataset file found in parent directory: %s', parent_dir);
        else
            % Load the first dataset file found (should only be one per project)
            dataset_file = fullfile(parent_dir, dataset_files(1).name);
            fprintf('    Loading dataset: %s\n', dataset_files(1).name);
            
            dataset_loaded = load(dataset_file);
            
            if ~isfield(dataset_loaded, 'data')
                warning('Dataset file does not contain data struct');
            else
                % Merge dataset data into current data
                % Dataset should have data.synGen
                % Run file has data.mcsbd_slice or data.mcsbd_block
                if isfield(dataset_loaded.data, 'synGen')
                    data.synGen = dataset_loaded.data.synGen;
                    fprintf('      Merged data.synGen from dataset\n');
                end
                
                % Merge params (dataset params are base, run params may override)
                % In practice, there should be no conflicts since generation is sequential
                dataset_params = dataset_loaded.params;
                param_fields = fieldnames(dataset_params);
                for i = 1:length(param_fields)
                    field = param_fields{i};
                    if ~isfield(params, field)
                        params.(field) = dataset_params.(field);
                    end
                end
                
                % Use dataset log as base, but keep run log as primary
                % Store dataset log in history
                load_history.dataset_log = dataset_loaded.log;
            end
            
            load_history.files_loaded{end+1} = dataset_file;
            load_history.load_order = [dataset_files(1).name, load_history.load_order];
        end
    end
    
    % Build meta struct
    meta = struct();
    meta.load_history = load_history;
    meta.loaded_file = workspace_file;
    meta.file_type = 'unknown';
    
    if is_run_file
        meta.file_type = 'run';
        meta.run_folder = current_folder;
        meta.run_file = file_name;
        meta.project_path = parent_dir;
    elseif is_dataset_file
        meta.file_type = 'dataset';
        meta.dataset_file = file_name;
        meta.project_path = file_dir;
    end
    
    % Display loaded data structure
    fprintf('\n');
    fprintf('  Loaded data structure:\n');
    data_fields = fieldnames(data);
    for i = 1:length(data_fields)
        fprintf('    data.%s\n', data_fields{i});
        if isstruct(data.(data_fields{i}))
            subfields = fieldnames(data.(data_fields{i}));
            for j = 1:min(5, length(subfields))
                fprintf('      .%s\n', subfields{j});
            end
            if length(subfields) > 5
                fprintf('      ... (%d more fields)\n', length(subfields) - 5);
            end
        end
    end
    
    fprintf('\n');
    fprintf('========================================================================\n');
    fprintf('  LOAD COMPLETE\n');
    fprintf('========================================================================\n');
    fprintf('  File type: %s\n', meta.file_type);
    fprintf('  Files loaded: %d\n', length(load_history.files_loaded));
    for i = 1:length(load_history.load_order)
        fprintf('    %d. %s\n', i, load_history.load_order{i});
    end
    fprintf('========================================================================\n');
end

