clc; clear;
run('../init_sbd');

%% Load pre-generated synthetic datasets
load('results\parallel results\synthetic_datasets_20241127_170302/synthetic_datasets_20241127_170302.mat');  
%%
% datasets = datasets(2);  % Choose subset if needed

%% Initialize kernels for all datasets
%A1_all = initialize_all_kernels(datasets);

%% Define system parameters and create configs
lambda1_range = [0.05, 0.2, 0.4];
mini_loop_range = [1, 3, 5];
maxIT = 30;

% Create parameter combinations
[lambda1_grid, mini_loop_grid] = meshgrid(lambda1_range, mini_loop_range);
param_combinations = [lambda1_grid(:), mini_loop_grid(:)];
num_param_combos = size(param_combinations, 1);

% Create parameter-specific configs upfront
root_dir = fileparts(pwd);
config_dir = fullfile(root_dir, 'config', 'worker_configs');  % Base config directory
if ~exist(config_dir, 'dir')
    mkdir(config_dir);
end

% Create config for each parameter combination
for m = 1:num_param_combos
    param_dir = fullfile(config_dir, sprintf('config_param_%d', m));
    if ~exist(param_dir, 'dir')
        mkdir(param_dir);
    end
    
    % Copy and update config files
    copyfile(fullfile(root_dir, 'config', 'Xsolve_config.mat'), ...
        fullfile(param_dir, 'Xsolve_config_tunable.mat'));
    copyfile(fullfile(root_dir, 'config', 'Asolve_config.mat'), ...
        fullfile(param_dir, 'Asolve_config_tunable.mat'));
    
    % Update configs with current parameters
    update_config(fullfile(param_dir, 'Xsolve_config_tunable.mat'), ...
        'MAXIT', param_combinations(m,2));
    update_config(fullfile(param_dir, 'Asolve_config_tunable.mat'), ...
        'options.maxiter', param_combinations(m,2));
end

%% Run parallel tests
% start initialize the parallel pool, check if already running
if isempty(gcp('nocreate'))
    parpool(9);
end

for n = 1:numel(datasets)   
    % print dataset index as indicator
    fprintf('Starting Dataset %d\n', n);
    % Broadcast parameters to workers
    % Extract dataset parameters that will be used in parfor
    dataset_Y = datasets(n).Y;
    dataset_kernel_size = datasets(n).params.kernel_size;
    dataset_X0 = datasets(n).X0;
    dataset_A0 = datasets(n).A0;
    dataset_A1 = A1_all{n};
    dataset_idx = n;

    % set parameters
    params.phase2 = false;
    params.X0 = dataset_X0;
    params.A0 = dataset_A0;
    params.Xsolve = 'FISTA';
    params.xpos = true;
    params.getbias = true;
    
    parfor param_idx = 1:num_param_combos
        % Use same path as creation
        param_dir = fullfile(config_dir, sprintf('config_param_%d', param_idx));
        try
            
            [Aout, Xout, bout, extras] = SBD_test_multi_parallel(...
                dataset_Y, ...
                dataset_kernel_size, ...
                params, ...
                dataset_A1, ...
                dataset_idx, ...
                param_combinations,...
                param_idx, ...
                maxIT);

            % Generate simpler filename based on dataset and parameter index
            filename = sprintf('SBD_parallel_dataset%d_param_idx%d.mat', ...
                n, param_idx);

            % Create a scalar structure (1x1) with all variables to save
            s = struct();  % Initialize empty scalar struct
            s(1).Aout = {Aout};  % Explicitly wrap in cell
            s(1).Xout = {Xout};
            s(1).bout = {bout};
            s(1).extras = {extras};
            s(1).param_combinations = param_combinations;
            s(1).param_idx = param_idx;
            s(1).dataset_A0 = {dataset_A0};
            s(1).dataset_X0 = {dataset_X0};

            % Save results immediately using -fromstruct option
            save(filename, '-fromstruct', s);
            
            fprintf('Results saved for dataset %d, param combo %d to: %s\n', ...
                n, param_idx, filename);

        catch ME
            % Enhanced error reporting
            fprintf('Error Details for dataset %d, param combo %d:\n', n, param_idx);
            fprintf('Error Message: %s\n', ME.message);
            fprintf('Error Stack:\n');
            for k = 1:length(ME.stack)
                fprintf('  File: %s\n  Line: %d\n  Function: %s\n', ...
                    ME.stack(k).file, ME.stack(k).line, ME.stack(k).name);
            end
        end
    end
end

%% Cleanup worker configs after completion
rmdir(config_dir, 's');