function [meta] = createProjectStructure(varargin)
%CREATEPROJECTSTRUCTURE Create project folder structure for synthetic data workflow
%
%   Creates hierarchical folder structure for organizing synthetic data
%   generation, algorithm runs, and visualization results.
%
%   [meta] = createProjectStructure(...)
%
%   INPUTS (Name-Value pairs):
%       'project_root'      - Root directory for projects (default: UI selection)
%                             If not provided, opens UI dialog to select directory
%       'timestamp'         - Timestamp string (default: auto-generated yyyymmdd_HHMMSS)
%       'create_folders'    - Actually create folders on disk (default: true)
%
%   OUTPUTS:
%       meta                - Metadata struct containing all path information:
%           .project_root   - Root directory for all projects
%           .project_name   - Project folder name (synthetic_<timestamp>)
%           .project_path   - Full path to project folder
%           .timestamp      - Timestamp used for project naming
%
%   FOLDER STRUCTURE CREATED (Nested Structure):
%       \MT-SBD-STM\projects\synthetic_<timestamp>\
%           ├─ auto\                        (dataset folder for auto initialization)
%           │  ├─ auto.mat                  (dataset file)
%           │  ├─ auto_LOGfile.txt          (dataset log file)
%           │  ├─ run01\                    (run folders)
%           │  │  ├─ run01.mat
%           │  │  ├─ run01_LOGfile.txt
%           │  │  └─ visualization01\      (visualization folders)
%           │  │     └─ (figure files)
%           │  ├─ run02\
%           │  │  └─ ...
%           │  └─ ...
%           ├─ manual01\                    (dataset folder for manual variant 1)
%           │  ├─ manual01.mat
%           │  ├─ manual01_LOGfile.txt
%           │  └─ run01\
%           │     └─ ...
%           └─ manual02\                    (dataset folder for manual variant 2)
%              ├─ manual02.mat
%              ├─ manual02_LOGfile.txt
%              └─ run01\
%                 └─ ...
%
%   DESCRIPTION:
%       This helper function creates the standardized nested folder structure for
%       the three-phase workflow (pre-run, post-run, visualization). It:
%       - Opens UI dialog to select project root directory (if not provided)
%       - Generates unique project folder with timestamp
%       - Creates project root and project folder on disk (if requested)
%       - Dataset folders are created by saveDataset.m when saving datasets
%       - Run folders are created by saveRun.m within dataset folders
%       - Visualization folders are created by saveVisualization.m within run folders
%       - Returns metadata struct with all paths for use by save functions
%
%   NOTE: This function only creates the project folder. Dataset, run, and
%   visualization folders are created automatically by their respective save
%   functions when needed.
%
%   EXAMPLE:
%       % Create project structure with UI selection (default)
%       meta = createProjectStructure();
%
%       % Create structure with specified project root
%       meta = createProjectStructure('project_root', 'C:/MyProjects');
%
%       % Create structure with custom timestamp
%       meta = createProjectStructure('timestamp', '20251031_120000');
%
%       % Create structure without making folders (for testing)
%       meta = createProjectStructure('create_folders', false);
%
%   See also: saveDataset, saveRun, loadSequential

    % Parse input arguments
    p = inputParser;
    addParameter(p, 'project_root', '', @ischar);
    addParameter(p, 'timestamp', '', @ischar);
    addParameter(p, 'create_folders', true, @islogical);
    parse(p, varargin{:});
    
    % Extract parameters
    project_root = p.Results.project_root;
    timestamp = p.Results.timestamp;
    create_folders = p.Results.create_folders;
    
    % Determine project root directory
    if isempty(project_root)
        % Use UI to select project root directory
        fprintf('Select project root directory...\n');
        default_path = pwd;
        project_root = uigetdir(default_path, 'Select Project Root Directory');
        
        if isequal(project_root, 0)
            error('No directory selected. Project creation cancelled.');
        end
        
        fprintf('  Selected project root: %s\n', project_root);
    end
    
    % Generate timestamp if not provided
    if isempty(timestamp)
        timestamp = datestr(now, 'yyyymmdd_HHMMSS');
    end
    
    % Create project folder name
    project_name = sprintf('synthetic_%s', timestamp);
    project_path = fullfile(project_root, project_name);
    
    % Create folders if requested
    if create_folders
        % Create project root if it doesn't exist
        if ~exist(project_root, 'dir')
            fprintf('  Creating projects root directory: %s\n', project_root);
            mkdir(project_root);
        end
        
        % Create project folder
        if ~exist(project_path, 'dir')
            fprintf('  Creating project folder: %s\n', project_path);
            mkdir(project_path);
        else
            warning('Project folder already exists: %s', project_path);
        end
    end
    
    % Build metadata struct
    meta = struct();
    meta.project_root = project_root;
    meta.project_name = project_name;
    meta.project_path = project_path;
    meta.timestamp = timestamp;
    
    fprintf('  Project structure initialized: %s\n', project_name);
end

