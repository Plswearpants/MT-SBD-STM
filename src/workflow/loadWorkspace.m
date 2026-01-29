function [log, data, params, meta] = loadWorkspace(varargin)
%LOADWORKSPACE Load saved workspace with memory overwrite protection
%
%   Loads a previously saved analysis state including data, parameters,
%   log configuration, and metadata. Uses a file browser UI to select a
%   .mat file and automatically finds the corresponding log file. Prompts
%   user to handle current workspace before loading (which OVERWRITES memory).
%
%   [log, data, params, meta] = loadWorkspace(...)
%
%   INPUTS (all optional Name-Value pairs):
%       'current_log'       - Current log struct (for memory protection)
%       'current_data'      - Current data struct (for save prompt)
%       'current_params'    - Current params struct (for save prompt)
%       'current_meta'      - Current meta struct (for context)
%
%       'workspace_file'    - Full path to .mat file to load (default: UI selection)
%                             If not provided, will open file browser dialog
%       'title'             - Title for file browser dialog (default: 'Select Workspace File')
%       'verify_log'        - Verify log file naming convention (default: true)
%       'track_selection'   - Store selected file in meta struct (default: false)
%
%   OUTPUTS:
%       log                 - Loaded log struct (REPLACES current_log)
%       data                - Loaded data struct (REPLACES current_data)
%       params              - Loaded params struct (REPLACES current_params)
%       meta                - Loaded/updated meta struct
%
%   DESCRIPTION:
%       This wrapper loads a complete analysis snapshot with safety checks:
%       - Simple file browser UI to select .mat file
%       - Automatically finds and verifies corresponding log file
%       - UI to select directory and name for NEW log file for current session
%       - Appends load entry to ORIGINAL log file (from loaded workspace)
%       - Updates log struct to point to new log file (so subsequent blocks log to new file)
%       - Memory overwrite protection (delete/rename/cancel existing workspace)
%       - Validates naming: <name>.mat and <name>_LOGfile.txt must exist in same directory
%       - Handles meta struct merging and selection tracking
%       - Loading ALWAYS overwrites log, data, params in memory!
%
%   EXAMPLE:
%       % First time load (file browser)
%       [log, data, params, meta] = loadWorkspace();
%
%       % Interactive load with existing workspace
%       [log, data, params, meta] = loadWorkspace('current_log', log, ...
%           'current_data', data, 'current_params', params, 'current_meta', meta);
%
%       % Load from specific file
%       [log, data, params, meta] = loadWorkspace('workspace_file', ...
%           'C:/results/experiment_001_20251030_143052.mat');
%
%       % Load with file browser (default)
%       [log, data, params, meta] = loadWorkspace('current_meta', meta, ...
%           'title', 'LOAD DATASET FOR ALGORITHM RUN', 'track_selection', true);
%
%   See also: saveDataset, saveRun, load, uigetfile

    % Parse input arguments
    p = inputParser;
    addParameter(p, 'current_log', struct(), @isstruct);
    addParameter(p, 'current_data', struct(), @isstruct);
    addParameter(p, 'current_params', struct(), @isstruct);
    addParameter(p, 'current_meta', struct(), @isstruct);
    addParameter(p, 'workspace_file', '', @ischar);
    addParameter(p, 'title', 'Select Workspace File', @ischar);
    addParameter(p, 'verify_log', true, @islogical);
    addParameter(p, 'track_selection', false, @islogical);
    parse(p, varargin{:});
    
    % Extract parameters
    current_log = p.Results.current_log;
    current_data = p.Results.current_data;
    current_params = p.Results.current_params;
    current_meta = p.Results.current_meta;
    workspace_file = p.Results.workspace_file;
    title_text = p.Results.title;
    verify_log = p.Results.verify_log;
    track_selection = p.Results.track_selection;
    
    % Determine workspace file - use UI if not provided
    if isempty(workspace_file)
        % Use file browser to select .mat file
        fprintf('Select workspace .mat file to load...\n');
        [filename, pathname] = uigetfile('*.mat', title_text, pwd);
        if isequal(filename, 0)
            error('No file selected. Load cancelled.');
        end
        workspace_file = fullfile(pathname, filename);
    end
    
    % Verify file exists
    if ~exist(workspace_file, 'file')
        error('Workspace file not found: %s', workspace_file);
    end
    
    % Extract workspace name and directory from file
    [workspace_dir, workspace_name, ~] = fileparts(workspace_file);
    
    % Find and verify matching log file exists
    log_file_path = fullfile(workspace_dir, [workspace_name '_LOGfile.txt']);
    if verify_log
        if ~exist(log_file_path, 'file')
            error(['Log file not found: %s\n' ...
                   'The log file must have the same name as the workspace:\n' ...
                   '  Expected: %s_LOGfile.txt\n' ...
                   '  Workspace: %s.mat\n' ...
                   'Please ensure naming convention is consistent.'], ...
                   log_file_path, workspace_name, workspace_name);
        else
            fprintf('Found corresponding log file: %s\n', log_file_path);
        end
    end
    
    % Check if current workspace exists and prompt user
    % Loading will OVERWRITE current log, data, params in memory!
    % Consider workspace exists if log has .file field (not empty struct)
    has_current_workspace = isstruct(current_log) && ...
                           isfield(current_log, 'file') && ...
                           ~isempty(current_log.file);
    
    if has_current_workspace
        fprintf('Current workspace detected in memory.\n')
        fprintf('Loading will overwrite current workspace in memory.\n');
        fprintf('WARNING: Current workspace is NOT automatically saved.\n');
        fprintf('If you need to save it, cancel this operation and use saveDataset or saveRun first.\n\n');
        
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
    
    % Load workspace
    fprintf('Loading workspace: %s\n', workspace_file);
    loaded = load(workspace_file);
    
    % Extract main structs (with fallback for different formats)
    if isfield(loaded, 'log')
        log = loaded.log;
    else
        error('Workspace file does not contain log struct');
    end
    
    if isfield(loaded, 'data')
        data = loaded.data;
    else
        error('Workspace file does not contain data struct');
    end
    
    if isfield(loaded, 'params')
        params = loaded.params;
    else
        error('Workspace file does not contain params struct');
    end
    
    % Handle meta struct
    if isfield(loaded, 'meta')
        % Merge loaded meta with current_meta (loaded takes precedence)
        meta_input_fields = fieldnames(current_meta);
        meta = loaded.meta;
        % Preserve fields from current_meta that aren't in loaded meta
        for i = 1:length(meta_input_fields)
            field = meta_input_fields{i};
            if ~isfield(meta, field)
                meta.(field) = current_meta.(field);
            end
        end
    else
        % No meta in file, use current_meta
        meta = current_meta;
    end
    
    % Track selection if requested
    if track_selection
        [~, filename_only, ~] = fileparts(workspace_file);
        meta.selected_dataset_file = [filename_only '.mat'];
        meta.selected_dataset_path = workspace_file;
    end
    
    % Store original log file info before creating new log file
    original_log_path = log.path;
    original_log_file = log.file;
    original_log_file_path = fullfile(original_log_path, [original_log_file '_LOGfile.txt']);
    
    % Verify original log file exists
    if verify_log
        if exist(original_log_file_path, 'file')
            fprintf('  ✓ Original log file verified: %s\n', original_log_file_path);
        else
            warning('Original log file not found: %s', original_log_file_path);
        end
    end
    
    % UI to select directory for new log file
    fprintf('\nSelect directory to save new log file...\n');
    % Suggest default path (original log path, or workspace directory, or current directory)
    default_path = original_log_path;
    if ~exist(default_path, 'dir')
        default_path = workspace_dir;
    end
    if ~exist(default_path, 'dir')
        default_path = pwd;
    end
    
    new_log_path = uigetdir(default_path, 'Select Directory for New Log File');
    if isequal(new_log_path, 0)
        error('No directory selected. Load cancelled.');
    end
    
    % Prompt user for log file name (use workspace name as default)
    default_log_name = sprintf('%s_loaded_%s', workspace_name, datestr(now, 'yyyymmdd_HHMMSS'));
    prompt = sprintf("Enter log file name (without _LOGfile.txt suffix) [%s]: ", default_log_name);
    log_file_input = input(prompt, 's');
    
    if isempty(log_file_input)
        log_file_input = default_log_name;
    end
    
    % Ensure log file name is char (not string object)
    log_file_input = char(log_file_input);
    
    % Ensure unique log file name
    new_log_file = log_file_input;
    new_log_file_path = fullfile(new_log_path, [new_log_file '_LOGfile.txt']);
    n = 1;
    while exist(new_log_file_path, 'file')
        new_log_file = sprintf('%s_%03d', log_file_input, n);
        new_log_file_path = fullfile(new_log_path, [new_log_file '_LOGfile.txt']);
        n = n + 1;
    end
    
    if n > 1
        fprintf('  Log file name already exists, using: %s\n', new_log_file);
    end
    
    % Initialize new log file with loading entry as first entry
    fprintf('  Creating new log file for current session: %s\n', new_log_file);
    % Initialize log file (creates header)
    logUsedBlocks(new_log_path, new_log_file, "LD01A", "Initializing log file", 1);
    % Append first entry: workspace loading (minimal: workspace path and new log path)
    LOGcomment = sprintf("Loaded from: %s | New log in: %s/%s_LOGfile.txt", workspace_file, new_log_path, new_log_file);
    LOGcomment = logUsedBlocks(new_log_path, new_log_file, "LD01A", LOGcomment, 0);
    
    % Update log struct to point to new log file
    log.path = new_log_path;
    log.file = new_log_file;
    
    % Display summary
    fprintf('Load complete.\n');
    fprintf('  Original log file: %s/%s_LOGfile.txt\n', original_log_path, original_log_file);
    fprintf('  New log file: %s/%s_LOGfile.txt\n', new_log_path, new_log_file);
    fprintf('\n');
    
end

