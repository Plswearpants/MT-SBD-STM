function [meta] = saveDataset(log, data, params, meta, varargin)
%SAVEDATASET Save synthetic data generation results (pre-run phase)
%
%   Saves the results of synthetic data generation and kernel initialization
%   with automatic naming based on initialization method (auto vs manual).
%
%   [meta] = saveDataset(log, data, params, meta, ...)
%
%   INPUTS:
%       log                 - Log struct with .path and .file fields
%       data                - Data struct with synGen field
%       params              - Parameter struct with all settings
%       meta                - Metadata struct from createProjectStructure
%
%       OPTIONAL (Name-Value pairs):
%       'compression'       - Use compression (-v7.3 for large files) (default: false)
%
%   OUTPUTS:
%       meta                - Updated metadata struct with:
%                             .dataset_file - Name of saved dataset file (e.g., 'auto', 'manual01')
%
%   SAVED FILES (Nested Structure):
%       <dataset>/<dataset>.mat              - MATLAB workspace with log, data, params
%       <dataset>/<dataset>_LOGfile.txt      - Copy of log file
%       Example: auto/auto.mat, manual01/manual01.mat
%
%   NAMING LOGIC:
%       - Auto initialization: 'auto.mat' (deterministic, no counter)
%       - Manual initialization: 'manual01.mat', 'manual02.mat', etc.
%         (sequential counter increments for each manual variant)
%
%   DESCRIPTION:
%       This wrapper saves the pre-run phase results (synthetic data + kernel
%       initialization) in nested structure. It:
%       - Determines dataset name from params.initialization.init_method
%       - Creates dataset folder: project_path/<dataset>/
%       - Saves to: project_path/<dataset>/<dataset>.mat
%       - Copies log file to dataset folder
%       - Updates meta.dataset_file (minimal - paths computed on-the-fly)
%
%   EXAMPLE:
%       % Save auto-initialized dataset
%       meta = saveDataset(log, data, params, meta);
%
%       % Save with compression
%       meta = saveDataset(log, data, params, meta, 'compression', true);
%
%   See also: createProjectStructure, saveRun, loadSequential

    % Parse input arguments
    p = inputParser;
    addRequired(p, 'log', @isstruct);
    addRequired(p, 'data', @isstruct);
    addRequired(p, 'params', @isstruct);
    addRequired(p, 'meta', @isstruct);
    addParameter(p, 'compression', false, @islogical);
    parse(p, log, data, params, meta, varargin{:});
    
    % Extract parameters
    compression = p.Results.compression;
    
    fprintf('Saving dataset (pre-run phase)...\n');
    
    % Validate inputs
    if ~isfield(data, 'synGen') || ~isstruct(data.synGen)
        error('Data struct must contain synGen field');
    end
    if ~isfield(meta, 'project_path')
        error('Meta struct must contain project_path field');
    end
    
    % Extract params from hierarchical structure if needed
    if isfield(params, 'synGen')
        params_temp = organizeParams(params, 'extract');
    else
        params_temp = params;
    end
    
    % Determine dataset type and file name from params.init_method
    if ~isfield(params_temp, 'init_method')
        error('Params struct must contain init_method field (from autoInitializeKernels or initializeKernelsRef)');
    end
    init_method = params_temp.init_method;
    
    if strcmp(init_method, 'auto')
        dataset_file = 'auto';
        fprintf('  Dataset type: auto\n');
    elseif startsWith(init_method, 'manual')
        dataset_file = init_method;  % 'manual01', 'manual02', etc.
        fprintf('  Dataset type: %s\n', init_method);
    else
        error('Unknown initialization method: %s', init_method);
    end
    
    % Create dataset folder (nested structure)
    dataset_folder = fullfile(meta.project_path, dataset_file);
    if ~exist(dataset_folder, 'dir')
        fprintf('  Creating dataset folder: %s\n', dataset_file);
        mkdir(dataset_folder);
    end
    
    % Build file paths (nested structure)
    mat_file = char(fullfile(dataset_folder, [dataset_file '.mat']));
    
    % Handle log.file (convert cell to char if necessary)
    if iscell(log.file)
        log_file_str = char(log.file{1});
    else
        log_file_str = char(log.file);
    end
    
    log_source = fullfile(log.path, [log_file_str '_LOGfile.txt']);
    log_dest = fullfile(dataset_folder, [dataset_file '_LOGfile.txt']);
    
    % Save workspace
    if compression
        fprintf('Saving with compression (v7.3)...\n');
        save(mat_file, 'log', 'data', 'params', '-v7.3');
    else
        fprintf('Saving without compression...\n');
        save(mat_file, 'log', 'data', 'params');
    end
    fprintf('Dataset saved: %s\n', [dataset_file '.mat']);
    
    % Copy log file
    if exist(log_source, 'file')
        copyfile(log_source, log_dest);
        fprintf('Log file copied: %s\n', [dataset_file '_LOGfile.txt']);
    else
        warning('Log file not found: %s', log_source);
    end
    
    % Update meta struct (minimal - only essential fields)
    meta.dataset_file = dataset_file;
    
    fprintf('Dataset save complete.\n');
    fprintf('  Location: %s/%s\n', meta.project_path, dataset_file);
end

