% Load and organize parallel test results for visualization
function [metrics] = load_parallel_results_old()
    % Initialize storage structure
    metrics = struct();
    
    % Get folder path from user
    default_path = fullfile(pwd, 'examples');
    folder_path = uigetdir(default_path, 'Select folder containing SBD parallel results');
    
    if folder_path == 0  % User canceled
        error('Folder selection canceled by user');
    end
    
    % Find and load synthetic dataset file first
    synthetic_files = dir(fullfile(folder_path, 'synthetic_datasets*.mat'));
    if ~isempty(synthetic_files)
        synthetic_data = load(fullfile(folder_path, synthetic_files(1).name));
        fprintf('Found synthetic dataset file: %s\n', synthetic_files(1).name);
        
        % Store dataset descriptions if available
        if isfield(synthetic_data, 'descriptions')
            metrics.dataset_descriptions = synthetic_data.descriptions;
        elseif isfield(synthetic_data.datasets, 'descriptions')
            metrics.dataset_descriptions = {synthetic_data.datasets.descriptions};
        end
    else
        warning('No synthetic dataset file found for descriptions');
    end
    
    % Get list of all result files
    files = dir(fullfile(folder_path, 'SBD_parallel_dataset*.mat'));
    
    if isempty(files)
        error('No SBD parallel result files found in selected folder');
    end
    
    % Display loading information
    fprintf('Found %d result files in:\n%s\n\n', length(files), folder_path);
    
    % Extract unique dataset numbers and parameter indices
    dataset_nums = [];
    param_indices = [];
    for i = 1:length(files)
        parts = split(files(i).name, {'_', '.'});
        dataset_nums = [dataset_nums; str2double(parts{3}(8:end))];
        param_indices = [param_indices; str2double(parts{5}(4:end))];
    end
    unique_datasets = unique(dataset_nums);
    num_params = max(param_indices);
    
    fprintf('Found data for %d unique datasets\n', length(unique_datasets));
    fprintf('Number of parameter combinations: %d\n\n', num_params);
    
    % Initialize arrays for metrics
    num_datasets = length(unique_datasets);
    metrics.kernel_quality_trajectory = cell(num_datasets, num_params);  % Store full trajectory
    metrics.activation_accuracy_trajectory = cell(num_datasets, num_params);  % Store full trajectory
    metrics.kernel_quality_final = nan(num_datasets, num_params);  % Store final values
    metrics.activation_accuracy_final = nan(num_datasets, num_params);  % Store final values
    metrics.runtime = nan(num_datasets, num_params);
    metrics.residuals = cell(num_datasets, num_params);
    metrics.relative_changes = cell(num_datasets, num_params);
    metrics.noise_levels = nan(num_datasets, 1);
    metrics.param_combinations = [];
    metrics.demixing_score = nan(num_datasets, num_params);
    metrics.demixing_matrices = cell(num_datasets, num_params);
    
    % Initialize progress bar
    fprintf('Loading and processing files:\n');
    progress_str = '';
    
    % Load each file and extract metrics
    for i = 1:length(files)
        % Update progress bar
        fprintf(repmat('\b', 1, length(progress_str)));
        progress_str = sprintf('Processing file %d/%d...', i, length(files));
        fprintf(progress_str);
        
        % Load result file
        data = load(fullfile(folder_path, files(i).name));
        
        % Extract dataset number and parameter index
        parts = split(files(i).name, {'_', '.'});
        dataset_idx = find(unique_datasets == str2double(parts{3}(8:end)));
        param_idx = str2double(parts{5}(4:end));
        
        % Load and check parameter combinations
        if isfield(data, 'param_combinations')
            current_params = data.param_combinations;
        elseif isfield(data, 's') && isfield(data.s, 'param_combinations')
            current_params = data.s.param_combinations;
        else
            error('Could not find param_combinations in data file');
        end
        
        % Store or append parameter combinations
        if isempty(metrics.param_combinations)
            metrics.param_combinations = current_params;
            fprintf('Initial parameter combinations loaded\n');
        else
            % Check if new combinations need to be appended
            new_combinations = [];
            for j = 1:size(current_params, 1)
                % Check if this combination already exists
                if ~any(all(abs(metrics.param_combinations - current_params(j,:)) < 1e-10, 2))
                    new_combinations = [new_combinations; current_params(j,:)];
                end
            end
            
            % Append new combinations if found
            if ~isempty(new_combinations)
                metrics.param_combinations = [metrics.param_combinations; new_combinations];
                fprintf('Added %d new parameter combinations\n', size(new_combinations, 1));
            end
        end
        
        % Extract metrics from extras
        extras = data.extras{1};  % Unwrap from cell
        
        % Store kernel quality (both trajectory and final value)
        if isfield(extras.phase1, 'kernel_quality_factors')
            kq = extras.phase1.kernel_quality_factors;
            metrics.kernel_quality_trajectory{dataset_idx, param_idx} = kq;  % Store full trajectory
            metrics.kernel_quality_final(dataset_idx, param_idx) = kq(end);  % Store final value
        end
        
        % Store activation accuracy (both trajectory and final value)
        if isfield(extras.phase1, 'activation_metrics')
            aa = extras.phase1.activation_metrics;
            metrics.activation_accuracy_trajectory{dataset_idx, param_idx} = aa;  % Store full trajectory
            metrics.activation_accuracy_final(dataset_idx, param_idx) = aa(end);  % Store final value
        end
        
        % Store runtime
        if isfield(extras, 'runtime')
            metrics.runtime(dataset_idx, param_idx) = extras.runtime;
        end
        
        % Calculate and store residuals
        if isfield(extras, 'residuals')
            metrics.residuals{dataset_idx, param_idx} = extras.residuals;
        end
        
        % Calculate relative changes
        if isfield(extras, 'residuals')
            res = extras.residuals;
            rel_changes = abs(diff(res)) ./ abs(res(1:end-1));
            metrics.relative_changes{dataset_idx, param_idx} = rel_changes;
        end
        
        % Store noise level (assuming its in extras)
        if isfield(extras, 'noise_level') && isnan(metrics.noise_levels(dataset_idx))
            metrics.noise_levels(dataset_idx) = extras.noise_level;
        end
        
        % Store demixing score and matrix
        if isfield(data, 'Xout')
            Xout = data.Xout{1};  % Assuming Xout is stored in a cell
            [demix_score, corr_matrix] = computeDemixingMetric(Xout);
            metrics.demixing_score(dataset_idx, param_idx) = demix_score;
            metrics.demixing_matrices{dataset_idx, param_idx} = corr_matrix;
        end
    end
    
    % Clear progress bar and print completion
    fprintf('\nData loading complete!\n\n');
    
    % Add metadata
    metrics.num_datasets = num_datasets;
    metrics.num_params = num_params;
    metrics.lambda1_values = unique(metrics.param_combinations(:,1));
    metrics.mini_loop_values = unique(metrics.param_combinations(:,2));
    
    % Print summary
    fprintf('Summary:\n');
    fprintf('- Number of datasets: %d\n', metrics.num_datasets);
    fprintf('- Lambda1 values: [%s]\n', join(string(metrics.lambda1_values), ', '));
    fprintf('- Mini-loop values: [%s]\n', join(string(metrics.mini_loop_values), ', '));

end