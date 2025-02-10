clc; clear;
run('../init_sbd');

%% Load pre-generated synthetic datasets
load('results\synthetic_datasets\synthetic_datasets_20250120_173448.mat');  

%% Initialize kernels for all datasets
fprintf('Initializing kernels for %d datasets...\n', length(datasets));
% Choose initialization method ('activation' or 'random')
init_method = 'random';

if strcmp(init_method, 'activation')
    A1_all = initialize_kernels_from_activations(datasets);
else
    A1_all = initialize_kernels_random(datasets);
end

%% Save initialized kernels to synthetic dataset
for i = 1:length(datasets)
    datasets(i).A1_init = A1_all{i};
end

% Save updated datasets with initialized kernels
[filepath, name, ext] = fileparts('results\synthetic_datasets\synthetic_datasets_20250120_173448.mat');
if strcmp(init_method, 'activation')
    init_type = '_init_activation';
else
    init_type = '_init_random';
end
new_filename = fullfile(filepath, [name, init_type, ext]);
save(new_filename, 'datasets');
fprintf('Saved initialized kernels to: %s\n', new_filename);

%% Visualize initialized kernels for random datasets
fprintf('Displaying initialized kernels for 4 random datasets...\n');
num_datasets = length(datasets);
num_to_display = min(4, num_datasets);
rand_dataset_indices = randperm(num_datasets, num_to_display);

figure('Name', 'Initialized Kernels for Random Datasets');
for idx = 1:num_to_display
    dataset_idx = rand_dataset_indices(idx);
    subplot(2, 2, idx);
    display_dataset_kernels(A1_all{dataset_idx}, datasets(dataset_idx), dataset_idx);
end
drawnow;

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
fprintf('Starting parallel processing of %d datasets\n', num_datasets);

% Pre-process all params_local structures outside parfor
params_all = cell(num_datasets, 1);
for n = 1:num_datasets
    params_all{n} = struct();
    params_all{n}.phase2 = false;
    params_all{n}.Xsolve = 'FISTA';
    params_all{n}.xpos = true;
    params_all{n}.getbias = true;
    params_all{n}.X0 = datasets(n).X0;
    params_all{n}.A0 = datasets(n).A0;
end

% Parallel loop over datasets
parfor n = 1:num_datasets
    try
        % Minimal progress indicator
        fprintf('Processing dataset %d/%d\n', n, num_datasets);
        
        % Extract dataset parameters
        dataset_Y = datasets(n).Y;
        dataset_kernel_size = datasets(n).params.kernel_size;
        dataset_A1 = A1_all{n};
        params_local = params_all{n};  % Get pre-processed params

        % Process dataset
        [Aout, Xout, bout, extras] = SBD_test_multi_parallel(...
            dataset_Y, ...
            dataset_kernel_size, ...
            params_local, ...
            dataset_A1, ...
            n, ...  % dataset_idx
            param_combinations,...
            param_idx, ...
            maxIT);

        % Save results with minimal console output
        filename = sprintf('SBD_parallel_dataset%d_optimal.mat', n);
        save_results(filename, Aout, Xout, bout, extras, param_combinations, ...
            param_idx, params_local.A0, params_local.X0);
        
    catch ME
        fprintf('Error in dataset %d: %s\n', n, ME.message);
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

function A1_all = initialize_kernels_random(datasets)
    % Initialize kernels using random Gaussian noise
    % Input:
    %   datasets: 1×N struct array with fields Y, X0, params
    % Output:
    %   A1_all: 1×N cell array of initial kernel guesses
    
    num_datasets = length(datasets);
    A1_all = cell(1, num_datasets);
    
    for i = 1:num_datasets
        kernel_sizes = datasets(i).params.kernel_size;  % [2×2] matrix of kernel sizes
        num_kernels = size(datasets(i).X0, 3);
        
        % Initialize kernels for this dataset
        A1 = cell(1, num_kernels);
        
        for k = 1:num_kernels
            % Get kernel size for this kernel
            kernel_height = kernel_sizes(k,1);
            kernel_width = kernel_sizes(k,2);
            
            % Generate random Gaussian noise kernel
            A1{k} = randn(kernel_height, kernel_width);
            
            % Normalize the kernel
            A1{k} = proj2oblique(A1{k});
            
            % Debug information
            fprintf('Dataset %d, Kernel %d: Random Gaussian kernel of size %dx%d\n', ...
                i, k, kernel_height, kernel_width);
        end
        
        A1_all{i} = A1;
    end
end

function display_dataset_kernels(A1_set, dataset, dataset_idx)
    % Display all kernels for a single dataset horizontally
    % Input:
    %   A1_set: cell array of initialized kernels for this dataset
    %   dataset: struct containing dataset information including kernel sizes
    %   dataset_idx: index of the dataset for display purposes
    
    kernel_sizes = dataset.params.kernel_size;
    num_kernels = length(A1_set);
    
    % Convert cell array to matrix for display
    max_height = max(kernel_sizes(:,1));
    max_width = max(kernel_sizes(:,2));
    kernel_matrix = zeros(max_height, num_kernels * max_width);
    
    % Arrange kernels horizontally with padding
    for k = 1:num_kernels
        current_kernel = A1_set{k};
        h_pad = floor((max_height - size(current_kernel,1))/2);
        w_pad = floor((max_width - size(current_kernel,2))/2);
        start_col = (k-1)*max_width + 1;
        kernel_matrix(h_pad+1:h_pad+size(current_kernel,1), ...
                    start_col+w_pad:start_col+w_pad+size(current_kernel,2)-1) = current_kernel;
    end
    
    % Display the set
    imagesc(kernel_matrix);
    colormap('gray');
    colorbar;
    title(sprintf('Dataset %d', dataset_idx));
    axis equal tight;
end

%% Helper function for saving results
function save_results(filename, Aout, Xout, bout, extras, param_combinations, ...
    param_idx, dataset_A0, dataset_X0)
    s = struct();
    s(1).Aout = {Aout};
    s(1).Xout = {Xout};
    s(1).bout = {bout};
    s(1).extras = {extras};
    s(1).param_combinations = param_combinations;
    s(1).param_idx = param_idx;
    s(1).dataset_A0 = {dataset_A0};
    s(1).dataset_X0 = {dataset_X0};
    save(filename, '-fromstruct', s);
end