function [metrics] = load_parallel_results()
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
    
    % Load first file to get parameter combinations
    first_data = load(fullfile(folder_path, files(1).name));
    if isfield(first_data, 'param_combinations')
        param_combinations = first_data.param_combinations;
    elseif isfield(first_data, 's') && isfield(first_data.s, 'param_combinations')
        param_combinations = first_data.s.param_combinations;
    else
        error('Could not find param_combinations in data file');
    end
    
    % Load and append any new parameter combinations from other files
    for i = 2:length(files)
        data = load(fullfile(folder_path, files(i).name));
        if isfield(data, 'param_combinations')
            current_params = data.param_combinations;
        elseif isfield(data, 's') && isfield(data.s, 'param_combinations')
            current_params = data.s.param_combinations;
        else
            continue;  % Skip if no param_combinations found
        end
        
        % Check for new combinations
        for j = 1:size(current_params, 1)
            % Check if this combination already exists
            if ~any(all(abs(param_combinations - current_params(j,:)) < 1e-10, 2))
                param_combinations = [param_combinations; current_params(j,:)];
                fprintf('Added new parameter combination: [%.3f, %d]\n', ...
                    current_params(j,1), current_params(j,2));
            end
        end
    end
    
    % Store final param_combinations in metrics
    metrics.param_combinations = param_combinations;
    
    % Get unique parameter values and dimensions
    lambda1_values = unique(param_combinations(:,1));
    mini_loop_values = unique(param_combinations(:,2));
    num_lambda1 = length(lambda1_values);
    num_mini_loop = length(mini_loop_values);
    num_datasets = length(unique_datasets);
    % print the number of datasets, lambda1 values, and mini-loop values
    fprintf('Number of datasets: %d\n', num_datasets);
    fprintf('Number of lambda1 values: %d\n', num_lambda1);
    fprintf('Number of mini-loop values: %d\n', num_mini_loop);

    fprintf('Found data for %d unique datasets\n', length(unique_datasets));
    fprintf('Number of parameter combinations: %d\n\n', num_lambda1 * num_mini_loop);
    
    % Initialize arrays for metrics with 3D structure: datasets × mini_loop × lambda1
    metrics.kernel_quality_trajectory = cell(num_datasets, num_mini_loop, num_lambda1);
    metrics.activation_accuracy_trajectory = cell(num_datasets, num_mini_loop, num_lambda1);
    metrics.kernel_quality_final = nan(num_datasets, num_mini_loop, num_lambda1);
    metrics.activation_accuracy_final = nan(num_datasets, num_mini_loop, num_lambda1);
    metrics.runtime = nan(num_datasets, num_mini_loop, num_lambda1);
    metrics.residuals = cell(num_datasets, num_mini_loop, num_lambda1);
    metrics.relative_changes = cell(num_datasets, num_mini_loop, num_lambda1);
    metrics.demixing_score = nan(num_datasets, num_mini_loop, num_lambda1);
    metrics.demixing_matrices = cell(num_datasets, num_mini_loop, num_lambda1);
    metrics.noise_levels = nan(num_datasets, 1);  % This remains 1D as its per dataset
    
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
        
        % Get current parameter values - using local param_combinations
        current_lambda1 = param_combinations(param_idx, 1);
        current_mini_loop = param_combinations(param_idx, 2);
        
        % Find indices in the 3D array
        lambda1_idx = find(lambda1_values == current_lambda1);
        mini_loop_idx = find(mini_loop_values == current_mini_loop);
        
        % Extract metrics from extras
        extras = data.extras{1};
        
        % Store metrics in 3D arrays
        if isfield(extras.phase1, 'kernel_quality_factors')
            kq = extras.phase1.kernel_quality_factors;
            metrics.kernel_quality_trajectory{dataset_idx, mini_loop_idx, lambda1_idx} = kq;
            metrics.kernel_quality_final(dataset_idx, mini_loop_idx, lambda1_idx) = prod(kq(end,:));
        end
        
        if isfield(extras.phase1, 'activation_metrics')
            aa = extras.phase1.activation_metrics;
            metrics.activation_accuracy_trajectory{dataset_idx, mini_loop_idx, lambda1_idx} = aa;
            metrics.activation_accuracy_final(dataset_idx, mini_loop_idx, lambda1_idx) = prod(aa(end,:));
        end
        
        if isfield(extras, 'runtime')
            metrics.runtime(dataset_idx, mini_loop_idx, lambda1_idx) = extras.runtime;
        end
        
        if isfield(extras, 'residuals')
            metrics.residuals{dataset_idx, mini_loop_idx, lambda1_idx} = extras.residuals;
            
            % Calculate relative changes
            res = extras.residuals;
            rel_changes = abs(diff(res)) ./ abs(res(1:end-1));
            metrics.relative_changes{dataset_idx, mini_loop_idx, lambda1_idx} = rel_changes;
        end
        
        if isfield(data, 'Xout')
            Xout = data.Xout{1};
            [demix_score, seperation_matrix] = computeDemixingMetric(Xout);
            metrics.demixing_score(dataset_idx, mini_loop_idx, lambda1_idx) = demix_score;
            metrics.demixing_matrices{dataset_idx, mini_loop_idx, lambda1_idx} = seperation_matrix;
        end
        
        % Store noise level (per dataset)
        if isfield(extras, 'noise_level') && isnan(metrics.noise_levels(dataset_idx))
            metrics.noise_levels(dataset_idx) = extras.noise_level;
        end
    end
    
    % Clear progress bar and print completion
    fprintf('\nData loading complete!\n\n');
    
    % Add metadata
    metrics.num_datasets = num_datasets;
    metrics.num_lambda1 = num_lambda1;
    metrics.num_mini_loop = num_mini_loop;
    metrics.lambda1_values = lambda1_values;
    metrics.mini_loop_values = mini_loop_values;
    
    % Print summary
    fprintf('Summary:\n');
    fprintf('- Number of datasets: %d\n', metrics.num_datasets);
    fprintf('- Lambda1 values: [%s]\n', join(string(metrics.lambda1_values), ', '));
    fprintf('- Mini-loop values: [%s]\n', join(string(metrics.mini_loop_values), ', '));
end