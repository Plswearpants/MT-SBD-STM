function [meta] = saveRun(log, data, params, meta, varargin)
%SAVERUN Save algorithm run results (post-run phase)
%
%   Saves the results of MT-SBD algorithm execution (single-slice or block)
%   in a runXX subfolder with sequential numbering.
%
%   [meta] = saveRun(log, data, params, meta, ...)
%
%   INPUTS:
%       log                 - Log struct with .path and .file fields
%       data                - Data struct with mcsbd_slice or mcsbd_block fields
%       params              - Parameter struct with all settings
%       meta                - Metadata struct with project_path
%
%       OPTIONAL (Name-Value pairs):
%       'compression'       - Use compression (-v7.3 for large files) (default: false)
%       'run_number'        - Specific run number to use (default: auto-increment)
%
%   OUTPUTS:
%       meta                - Updated metadata struct with:
%                             .run_number - Run number used
%                             .run_folder - Name of run folder (runXX)
%                             .run_path - Full path to run folder
%                             .run_file - Name of run file (runXX.mat)
%                             .run_file_path - Full path to run file
%                             .run_type - Type of run ('slice', 'block', or 'both')
%                             .selected_dataset_file - Dataset file used (e.g., 'auto.mat')
%
%   SAVED FILES:
%       \run_<dataset>_XX\run_<dataset>_XX.mat      - MATLAB workspace
%       \run_<dataset>_XX\run_<dataset>_XX_LOGfile.txt - Copy of log file
%
%   FOLDER STRUCTURE:
%       \MT-SBD-STM\projects\synthetic_<timestamp>\
%           ├─ auto.mat or manualXX.mat (dataset files)
%           ├─ run_auto_01\              (runs based on auto.mat)
%           │  ├─ run_auto_01.mat
%           │  └─ run_auto_01_LOGfile.txt
%           ├─ run_auto_02\
%           │  ├─ run_auto_02.mat
%           │  └─ run_auto_02_LOGfile.txt
%           ├─ run_manual01_01\          (runs based on manual01.mat)
%           │  ├─ run_manual01_01.mat
%           │  └─ run_manual01_01_LOGfile.txt
%           └─ run_manual02_01\          (runs based on manual02.mat)
%              ├─ run_manual02_01.mat
%              └─ run_manual02_01_LOGfile.txt
%
%   DESCRIPTION:
%       This wrapper saves the post-run phase results (algorithm execution).
%       It:
%       - Creates run_<dataset>_XX subfolder in project directory
%       - Auto-increments run number based on existing runs for that dataset
%       - Saves data, params, log, and meta structs to .mat file
%       - Copies log file to run folder
%       - Updates meta struct with run information
%       - Uses dataset name in folder/file names for easy identification
%
%   EXAMPLE:
%       % Save run with auto-increment
%       meta = saveRun(log, data, params, meta);
%
%       % Save specific run number with compression
%       meta = saveRun(log, data, params, meta, 'run_number', 5, 'compression', true);
%
%   See also: createProjectStructure, saveDataset, loadSequential

    % Parse input arguments
    p = inputParser;
    addRequired(p, 'log', @isstruct);
    addRequired(p, 'data', @isstruct);
    addRequired(p, 'params', @isstruct);
    addRequired(p, 'meta', @isstruct);
    addParameter(p, 'compression', false, @islogical);
    addParameter(p, 'run_number', [], @(x) isempty(x) || isnumeric(x));
    parse(p, log, data, params, meta, varargin{:});
    
    % Extract parameters
    compression = p.Results.compression;
    run_number = p.Results.run_number;
    
    fprintf('Saving run results (post-run phase)...\n');
    
    % Validate inputs
    if ~isfield(meta, 'project_path')
        error('Meta struct must contain project_path field');
    end
    
    % Check that data contains algorithm results
    has_slice_results = isfield(data, 'mcsbd_slice') && isstruct(data.mcsbd_slice);
    has_block_results = isfield(data, 'mcsbd_block') && isstruct(data.mcsbd_block);
    
    if ~has_slice_results && ~has_block_results
        error('Data struct must contain mcsbd_slice or mcsbd_block field');
    end
    
    % Determine run type
    if has_slice_results && has_block_results
        run_type = 'both';
    elseif has_slice_results
        run_type = 'slice';
    else
        run_type = 'block';
    end
    fprintf('  Run type: %s\n', run_type);
    
    % Check which dataset was used (if tracked) and extract dataset name
    dataset_name = '';
    if isfield(meta, 'selected_dataset_file')
        fprintf('  Based on dataset: %s\n', meta.selected_dataset_file);
        % Extract dataset name without extension (e.g., 'auto', 'manual01')
        [~, dataset_name, ~] = fileparts(meta.selected_dataset_file);
    end
    
    % Determine run number based on dataset-specific runs
    if isempty(run_number)
        % Auto-increment: find highest existing run number for this dataset
        run_number = 1;
        
        if ~isempty(dataset_name)
            % Look for runs with this dataset name pattern
            run_pattern = sprintf('run_%s*', dataset_name);
            run_folders = dir(fullfile(meta.project_path, run_pattern));
            run_folders = run_folders([run_folders.isdir]);
            
            for i = 1:length(run_folders)
                folder_name = run_folders(i).name;
                % Extract number from 'run_<dataset>_XX' format
                expected_prefix = sprintf('run_%s_', dataset_name);
                if startsWith(folder_name, expected_prefix)
                    num_str = folder_name(length(expected_prefix)+1:end);
                    num = str2double(num_str);
                    if ~isnan(num) && num >= run_number
                        run_number = num + 1;
                    end
                end
            end
            fprintf('  Auto-incremented run number for %s: %02d\n', dataset_name, run_number);
        else
            % Fallback to generic run numbering if no dataset name
            run_folders = dir(fullfile(meta.project_path, 'run*'));
            run_folders = run_folders([run_folders.isdir]);
            
            for i = 1:length(run_folders)
                folder_name = run_folders(i).name;
                if startsWith(folder_name, 'run')
                    num_str = folder_name(4:end);
                    num = str2double(num_str);
                    if ~isnan(num) && num >= run_number
                        run_number = num + 1;
                    end
                end
            end
            fprintf('  Auto-incremented run number: %02d\n', run_number);
        end
    else
        fprintf('  Using specified run number: %02d\n', run_number);
    end
    
    % Create run folder name and path
    if ~isempty(dataset_name)
        run_folder = sprintf('run_%s_%02d', dataset_name, run_number);
    else
        run_folder = sprintf('run%02d', run_number);
    end
    run_path = fullfile(meta.project_path, run_folder);
    
    % Create run folder if it doesn't exist
    if ~exist(run_path, 'dir')
        fprintf('  Creating run folder: %s\n', run_folder);
        mkdir(run_path);
    else
        warning('Run folder already exists: %s (will overwrite)', run_folder);
    end
    
    % Build file paths (use same name as folder for consistency)
    run_file = run_folder;  % e.g., 'run_auto_01' or 'run_manual01_01'
    mat_file = char(fullfile(run_path, [run_file '.mat']));
    
    % Handle log.file (convert cell to char if necessary)
    if iscell(log.file)
        log_file_str = char(log.file{1});
    else
        log_file_str = char(log.file);
    end
    
    log_source = fullfile(log.path, [log_file_str '_LOGfile.txt']);
    log_dest = fullfile(run_path, [run_file '_LOGfile.txt']);
    
    % Save workspace (including meta struct)
    fprintf('  Saving run results...\n');
    if compression
        fprintf('    Saving with compression (v7.3)...\n');
        save(mat_file, 'log', 'data', 'params', 'meta', '-v7.3');
    else
        fprintf('    Saving...\n');
        save(mat_file, 'log', 'data', 'params', 'meta');
    end
    fprintf('    Run saved: %s\n', [run_file '.mat']);
    
    % Copy log file
    if exist(log_source, 'file')
        copyfile(log_source, log_dest);
        fprintf('    Log file copied: %s\n', [run_file '_LOGfile.txt']);
    else
        warning('Log file not found: %s', log_source);
    end
    
    % Update meta struct
    meta.run_number = run_number;
    meta.run_folder = run_folder;
    meta.run_path = run_path;
    meta.run_file = run_file;
    meta.run_file_path = mat_file;
    meta.run_type = run_type;
    
    fprintf('  Run save complete.\n');
    fprintf('  Location: %s\n', run_path);
end

