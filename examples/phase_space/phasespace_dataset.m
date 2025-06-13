clc; clear;
run('../init_sbd');

%% Load pre-generated synthetic datasets
load('results\synthetic_datasets\synthetic_datasets_20250611_182057.mat');  

%% Visualize pre-initialized kernels for random datasets
fprintf('Displaying pre-initialized kernels for 4 random datasets...\n');
num_datasets = length(datasets);
num_to_display = min(4, num_datasets);
rand_dataset_indices = randperm(num_datasets, num_to_display);

figure('Name', 'Pre-initialized Kernels for Random Datasets');
for idx = 1:num_to_display
    dataset_idx = rand_dataset_indices(idx);
    subplot(2, 2, idx);
    display_dataset_kernels(datasets(dataset_idx).A1, datasets(dataset_idx), dataset_idx);
end
drawnow;

%% Define optimal parameters
lambda1_optimal = 0.2;  % optimal regularization strength
mini_loop_optimal = 2;  % optimal iteration count
maxIT = 10;

% Create single parameter combination
param_combinations = [lambda1_optimal, mini_loop_optimal];
param_idx = 1;  

% Create config directory
root_dir = 'C:\Users\CAD\Documents\GitHub\MT-SBD-STM';
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
    params_all{n}.phase2 = true;  % Enable Phase 2
    params_all{n}.Xsolve = 'FISTA';
    params_all{n}.xpos = true;
    params_all{n}.getbias = true;
    params_all{n}.X0 = datasets(n).X0;
    params_all{n}.A0 = datasets(n).A0;
    
    % Add Phase 1 parameters
    params_all{n}.lambda1 = lambda1_optimal;  % Regularization parameter for Phase I
    
    % Add Phase 2 parameters
    params_all{n}.kplus = ceil(0.3 * datasets(n).params.kernel_size);  % Border padding for sphere lifting
    params_all{n}.lambda2 = params_all{n}.lambda1 * 0.5;  % FINAL reg. param. value for Phase II 
    params_all{n}.nrefine = 8;  % Number of refinements for Phase II
end

% Parallel loop over datasets
parfor n = 1:num_datasets
    try
        % Minimal progress indicator
        fprintf('Processing dataset %d/%d\n', n, num_datasets);
        
        % Extract dataset parameters
        dataset_Y = datasets(n).Y;
        dataset_kernel_size = datasets(n).params.kernel_size;
        dataset_A1 = datasets(n).A1;  % Use pre-initialized kernels
        params_local = params_all{n};  % Get pre-processed params

        % Process dataset
        [Aout, Xout, bout, extras] = SBD_test_multi_parallel(...
            dataset_Y, ...
            dataset_kernel_size, ...
            params_all{n}, ...
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