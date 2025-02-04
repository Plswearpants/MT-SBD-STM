clc; clear;
run('../init_sbd');

%% Load pre-generated synthetic datasets
load('results\parallel results\synthetic_datasets_20250114_101330/synthetic_datasets_20250114_101330.mat');  

%% Initialize kernels for all datasets
fprintf('Initializing kernels for %d datasets...\n', length(datasets));
A1_all = initialize_kernels_from_activations(datasets);

%% Define optimal parameters (no longer need parameter loops)
lambda1_optimal = 0.2;  % example value
mini_loop_optimal = 3;  % example value
maxIT = 30;

% Create single parameter combination
param_combinations = [lambda1_optimal, mini_loop_optimal];
param_idx = 1;  % Since we only have one combination now

% Create config directory
root_dir = fileparts(pwd);
config_dir = fullfile(root_dir, 'config', 'worker_configs');
if ~exist(config_dir, 'dir')
    mkdir(config_dir);
end

% Create single config with optimal parameters
param_dir = fullfile(config_dir, 'config_param_1');
if ~exist(param_dir, 'dir')
    mkdir(param_dir);
end

% Copy and update config files
copyfile(fullfile(root_dir, 'config', 'Xsolve_config.mat'), ...
    fullfile(param_dir, 'Xsolve_config_tunable.mat'));
copyfile(fullfile(root_dir, 'config', 'Asolve_config.mat'), ...
    fullfile(param_dir, 'Asolve_config_tunable.mat'));

% Update configs with optimal parameters
update_config(fullfile(param_dir, 'Xsolve_config_tunable.mat'), ...
    'MAXIT', mini_loop_optimal);
update_config(fullfile(param_dir, 'Asolve_config_tunable.mat'), ...
    'options.maxiter', mini_loop_optimal);

%% Initialize parallel pool
if isempty(gcp('nocreate'))
    parpool(9);  % Creates pool with 9 workers
end

%% Run parallel processing on datasets
num_datasets = numel(datasets);
fprintf('Processing %d datasets with %d workers\n', num_datasets, gcp().NumWorkers);

% Define params structure outside parfor loop
params = struct();
params.phase2 = false;
params.Xsolve = 'FISTA';
params.xpos = true;
params.getbias = true;

% Pre-allocate cell array for results
results = cell(num_datasets, 1);

% Parallel loop over datasets
parfor n = 1:num_datasets
    try
        % Print progress
        fprintf('Starting Dataset %d\n', n);
        
        % Extract dataset parameters
        dataset_Y = datasets(n).Y;
        dataset_kernel_size = datasets(n).params.kernel_size;
        dataset_X0 = datasets(n).X0;
        dataset_A0 = datasets(n).A0;
        dataset_A1 = A1_all{n};
        dataset_idx = n;

        % Add debug print to verify kernel sizes
        fprintf('Dataset %d kernel sizes:\n', n);
        for k = 1:size(dataset_kernel_size, 1)
            fprintf('  Kernel %d: [%d × %d]\n', k, dataset_kernel_size(k,1), dataset_kernel_size(k,2));
        end
        
        % Verify data format
        fprintf('Dataset %d - X0 size: [%s], num A0 kernels: %d\n', ...
            n, mat2str(size(dataset_X0)), length(dataset_A0));
        
        % Create local copy of params and add dataset-specific fields
        params_local = params;
        params_local.X0 = dataset_X0;  % Should be [height × width × num_kernels]
        params_local.A0 = dataset_A0;  % Should be cell array of kernels

        % Process dataset
        [Aout, Xout, bout, extras] = SBD_test_multi_parallel(...
            dataset_Y, ...
            dataset_kernel_size, ...
            params_local, ...
            dataset_A1, ...
            dataset_idx, ...
            param_combinations,...
            param_idx, ...
            maxIT);

        % Generate filename
        filename = sprintf('SBD_parallel_dataset%d_optimal.mat', n);

        % Create results structure
        s = struct();
        s(1).Aout = {Aout};
        s(1).Xout = {Xout};
        s(1).bout = {bout};
        s(1).extras = {extras};
        s(1).param_combinations = param_combinations;
        s(1).param_idx = param_idx;
        s(1).dataset_A0 = {dataset_A0};
        s(1).dataset_X0 = {dataset_X0};

        % Save results
        save(filename, '-fromstruct', s);
        
        fprintf('Results saved for dataset %d to: %s\n', n, filename);

    catch ME
        % Error reporting
        fprintf('Error processing dataset %d:\n', n);
        fprintf('Error Message: %s\n', ME.message);
        fprintf('Error Stack:\n');
        for k = 1:length(ME.stack)
            fprintf('  File: %s\n  Line: %d\n  Function: %s\n', ...
                ME.stack(k).file, ME.stack(k).line, ME.stack(k).name);
        end
    end
end

%% Cleanup worker configs
rmdir(config_dir, 's');

fprintf('Processing complete!\n');

%% Helper Functions
function A1_all = initialize_kernels_from_activations(datasets)
    % Initialize kernels using activation point closest to center for each dataset
    % Input:
    %   datasets: 1×N struct array with fields Y, X0, params
    % Output:
    %   A1_all: 1×N cell array of initial kernel guesses
    
    num_datasets = length(datasets);
    A1_all = cell(1, num_datasets);
    
    for i = 1:num_datasets
        % Get data for this dataset
        Y = datasets(i).Y;
        X0 = datasets(i).X0;
        kernel_sizes = datasets(i).params.kernel_size;  % [2×2] matrix of kernel sizes
        num_kernels = size(X0, 3);
        
        % Calculate image center
        center = round(size(Y)/2);
        
        % Initialize kernels for this dataset
        A1 = cell(1, num_kernels);
        
        for k = 1:num_kernels
            % Get kernel size for this kernel
            kernel_height = kernel_sizes(k,1);
            kernel_width = kernel_sizes(k,2);
            half_height = floor(kernel_height/2);
            half_width = floor(kernel_width/2);
            
            % Find all activation points
            [rows, cols] = find(X0(:,:,k));
            
            if ~isempty(rows)
                % Calculate distance to center for each activation point
                distances = sqrt((rows - center(1)).^2 + (cols - center(2)).^2);
                [~, idx] = min(distances);
                
                % Get closest point to center
                row = rows(idx);
                col = cols(idx);
                
                % Calculate initial window bounds
                row_start = row - half_height;
                row_end = row + half_height;
                col_start = col - half_width;
                col_end = col + half_width;
                
                % Adjust window to fit within image bounds
                if row_start < 1
                    shift = 1 - row_start;
                    row_start = 1;
                    row_end = row_end + shift;
                elseif row_end > size(Y,1)
                    shift = row_end - size(Y,1);
                    row_end = size(Y,1);
                    row_start = row_start - shift;
                end
                
                if col_start < 1
                    shift = 1 - col_start;
                    col_start = 1;
                    col_end = col_end + shift;
                elseif col_end > size(Y,2)
                    shift = col_end - size(Y,2);
                    col_end = size(Y,2);
                    col_start = col_start - shift;
                end
                
                % Extract and normalize initial kernel
                A1{k} = Y(row_start:row_end, col_start:col_end);
                
                % Resize to match exact kernel size if needed
                if size(A1{k},1) ~= kernel_height || size(A1{k},2) ~= kernel_width
                    A1{k} = imresize(A1{k}, [kernel_height, kernel_width]);
                end
                
                A1{k} = proj2oblique(A1{k});  % Normalize
                
                % Debug information
                fprintf('Dataset %d, Kernel %d: Window [%d:%d, %d:%d], Size: %dx%d (Target: %dx%d)\n', ...
                    i, k, row_start, row_end, col_start, col_end, ...
                    size(A1{k},1), size(A1{k},2), kernel_height, kernel_width);
            else
                warning('No activation found for kernel %d in dataset %d', k, i);
                % Fallback: use center of image
                A1{k} = Y(center(1)-half_height:center(1)+half_height, ...
                         center(2)-half_width:center(2)+half_width);
                
                % Resize to match exact kernel size if needed
                if size(A1{k},1) ~= kernel_height || size(A1{k},2) ~= kernel_width
                    A1{k} = imresize(A1{k}, [kernel_height, kernel_width]);
                end
                
                A1{k} = proj2oblique(A1{k});
            end
        end
        
        A1_all{i} = A1;
    end
end