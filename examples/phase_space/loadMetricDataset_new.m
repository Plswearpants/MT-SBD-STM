function dataset_metrics = loadMetricDataset_new()
    % Initialize storage structure
    dataset_metrics = struct();
    
    % Get folder path from user with error handling
    default_path = fullfile(pwd, 'examples');
    folder_path = uigetdir(default_path, 'Select folder containing SBD results');
    if folder_path == 0
        error('Folder selection canceled by user');
    end
    
    % First load synthetic dataset file
    synthetic_files = dir(fullfile(folder_path, 'synthetic_datasets*.mat'));
    if isempty(synthetic_files)
        error('No synthetic dataset file found');
    end

    % Load parameter information
    try
        synthetic_data = load(fullfile(folder_path, synthetic_files(1).name), 'datasets', 'param_sets', 'descriptions');
        param_sets = synthetic_data.param_sets;  % [SNR, theta_cap, Nobs]
    catch ME
        error('Failed to load synthetic dataset file: %s', ME.message);
    end

    % Create parameter mapping (hash table approach)
    [unique_snr, ~, snr_idx] = unique(param_sets(:,1));
    [unique_theta, ~, theta_idx] = unique(param_sets(:,2));
    [unique_nobs, ~, nobs_idx] = unique(param_sets(:,3));
    
    % Store parameter values
    dataset_metrics.SNR_values = unique_snr;
    dataset_metrics.theta_cap_values = unique_theta;
    dataset_metrics.Nobs_values = unique_nobs;
    
    % Create parameter lookup table
    param_lookup = zeros(size(param_sets, 1), 3);
    for i = 1:size(param_sets, 1)
        param_lookup(i,:) = [snr_idx(i), theta_idx(i), nobs_idx(i)];
    end
    
    % Initialize metric arrays with pre-allocated memory
    dims = [length(unique_snr), length(unique_theta), length(unique_nobs)];
    dataset_metrics.kernel_quality_trajectory = cell(dims);
    dataset_metrics.activation_similarity_trajectory = cell(dims);
    dataset_metrics.kernel_quality_final = nan(dims);
    dataset_metrics.activation_similarity_final = nan(dims);
    dataset_metrics.runtime = nan(dims);
    dataset_metrics.residuals = cell(dims);
    dataset_metrics.relative_changes = cell(dims);
    dataset_metrics.combined_activationScore = nan(dims);
    dataset_metrics.demixing_score = nan(dims);
    
    % Add storage for reconstruction data
    dataset_metrics.Y = cell(dims);          % Original observations from synthetic data
    dataset_metrics.Y_clean = cell(dims);    % Clean observations from synthetic data
    dataset_metrics.A0 = cell(dims);         % Ground truth kernels
    dataset_metrics.A0_noiseless = cell(dims); % noiseless ground truth kernels 
    dataset_metrics.A0_noise_normalized = cell(dims); % the best ground truth kernels the algorithm can get  
    dataset_metrics.X0 = cell(dims);         % Ground truth activations
    dataset_metrics.Aout = cell(dims);       % Reconstructed kernels
    dataset_metrics.Xout = cell(dims);       % Reconstructed activations
    dataset_metrics.bout = cell(dims);       % Bias terms
    dataset_metrics.extras = cell(dims);    % Extra metrics

    % Get list of result files
    result_files = dir(fullfile(folder_path, 'SBD_parallel_dataset*_optimal.mat'));
    if isempty(result_files)
        error('No SBD result files found');
    end
    
    % Process result files in batches for memory efficiency
    batch_size = 10;
    num_batches = ceil(length(result_files) / batch_size);
    
    for batch = 1:num_batches
        start_idx = (batch-1)*batch_size + 1;
        end_idx = min(batch*batch_size, length(result_files));
        
        % Process files in current batch
        for i = start_idx:end_idx
            % Extract dataset number efficiently
            [~, name] = fileparts(result_files(i).name);
            dataset_num = str2double(regexp(name, '\d+', 'match'));
            
            % Get pre-computed parameter indices
            indices = param_lookup(dataset_num,:);
            
            % Process current dataset
            try
                % Load result data
                data = load(fullfile(folder_path, result_files(i).name));
                
                % Store original and clean observations
                dataset_metrics.Y{indices(1), indices(2), indices(3)} = synthetic_data.datasets(dataset_num).Y;
                dataset_metrics.Y_clean{indices(1), indices(2), indices(3)} = synthetic_data.datasets(dataset_num).Y_clean;
                dataset_metrics.A0_noiseless{indices(1), indices(2), indices(3)} = synthetic_data.datasets(dataset_num).A0_noiseless;
                % Store ground truth data
                if isfield(data, 'dataset_A0')
                    dataset_metrics.A0{indices(1), indices(2), indices(3)} = data.dataset_A0{1};
                end
                if isfield(data, 'dataset_X0')
                    dataset_metrics.X0{indices(1), indices(2), indices(3)} = data.dataset_X0{1};
                end
                
                % Store reconstruction data directly without flipping check
                if isfield(data, 'Xout')
                    dataset_metrics.Xout{indices(1), indices(2), indices(3)} = data.Xout{1};
                end
                if isfield(data, 'Aout')
                    dataset_metrics.Aout{indices(1), indices(2), indices(3)} = data.Aout{1};
                end
                if isfield(data, 'bout')
                    dataset_metrics.bout{indices(1), indices(2), indices(3)} = data.bout{1};
                end
                
                % Calculate metrics for original order
                if isfield(data, 'Xout') && isfield(data, 'dataset_X0') && ...
                   ~isempty(data.Xout{1}) && ~isempty(data.dataset_X0{1}) && ...
                   isfield(data, 'Aout') && isfield(data, 'dataset_A0')
                    % Get quality metrics without flipping
                    %[~, quality_scores] = detect_kernel_flip(data.dataset_X0, data.Xout, data.dataset_A0, data.Aout, true);
                    kk = synthetic_data.datasets(dataset_num).A0_noiseless;
                    [~, quality_scores] = detect_kernel_flip(data.dataset_X0, data.Xout, kk, data.Aout, true);
                    % Store metrics from original order
                    kernel_quality_final = mean(quality_scores.no_flip.kernel_similarity);
                    dataset_metrics.kernel_quality_final(indices(1), indices(2), indices(3)) = kernel_quality_final;

                    act_sim_final = mean(quality_scores.no_flip.activation_similarity);
                    dataset_metrics.activation_similarity_final(indices(1), indices(2), indices(3)) = act_sim_final;

                    % Calculate demixing score
                    [demixing_score, ~] = computeDemixingMetric(data.Xout{1});
                    dataset_metrics.demixing_score(indices(1), indices(2), indices(3)) = demixing_score;
                    
                    % Calculate combined score
                    dataset_metrics.combined_activationScore(indices(1), indices(2), indices(3)) = ...
                        computeCombined_activationScore(demixing_score, act_sim_final);
                end
                
                % Load metrics from trajectories
                extras = data.extras{1};
                dataset_metrics.extras{indices(1), indices(2), indices(3)} = extras;
                if isstruct(extras) && isfield(extras, 'phase1')
                    % Store trajectories
                    if isfield(extras.phase1, 'kernel_quality_factors')
                        kq_traj = extras.phase1.kernel_quality_factors;
                        dataset_metrics.kernel_quality_trajectory{indices(1), indices(2), indices(3)} = kq_traj;
                    end
                    
                    if isfield(extras.phase1, 'activation_metrics')
                        act_traj = extras.phase1.activation_metrics;
                        dataset_metrics.activation_similarity_trajectory{indices(1), indices(2), indices(3)} = act_traj;
                    end
                    
                    % Store other phase1 metrics
                    if isfield(extras.phase1, 'residuals')
                        dataset_metrics.residuals{indices(1), indices(2), indices(3)} = extras.phase1.residuals;
                    end
                    if isfield(extras.phase1, 'relative_changes')
                        dataset_metrics.relative_changes{indices(1), indices(2), indices(3)} = extras.phase1.relative_changes;
                    end
                end
                
                % Store runtime
                if isfield(extras, 'runtime')
                    dataset_metrics.runtime(indices(1), indices(2), indices(3)) = extras.runtime;
                end
                
            catch ME
                warning('Error processing file %s: %s', result_files(i).name, ME.message);
                continue;
            end
        end
        
        % Display progress
        fprintf('Processed batch %d/%d\n', batch, num_batches);
    end
    
    % Print summary
    print_metrics_summary(dataset_metrics);
end

function print_metrics_summary(metrics)
    fprintf('\nDataset Parameters Summary:\n');
    fprintf('- SNR: %d values [%.2e to %.2e]\n', ...
        length(metrics.SNR_values), min(metrics.SNR_values), max(metrics.SNR_values));
    fprintf('- Theta cap: %d values [%.2e to %.2e]\n', ...
        length(metrics.theta_cap_values), min(metrics.theta_cap_values), max(metrics.theta_cap_values));
    fprintf('- Nobs: %d values [%.2f to %.2f]\n', ...
        length(metrics.Nobs_values), min(metrics.Nobs_values), max(metrics.Nobs_values));
    
    % Add summary of stored reconstruction data
    fprintf('\nReconstruction Data Summary:\n');
    non_empty_reconstructions = nnz(~cellfun(@isempty, metrics.Aout));
    total_points = numel(metrics.Aout);
    fprintf('- Stored reconstructions: %d/%d points (%.1f%%)\n', ...
        non_empty_reconstructions, total_points, ...
        100 * non_empty_reconstructions / total_points);
end 