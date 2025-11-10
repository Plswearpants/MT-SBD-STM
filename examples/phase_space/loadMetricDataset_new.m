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
        % Note: param_sets may have fewer rows than datasets if repetitions were used
        % We'll extract unique parameter values from the actual datasets
    catch ME
        error('Failed to load synthetic dataset file: %s', ME.message);
    end

    % Extract unique parameter values from datasets (handles repetitions)
    % Get all parameter values from datasets
    all_snr = [];
    all_theta = [];
    all_nobs = [];
    all_reps = [];
    for i = 1:length(synthetic_data.datasets)
        if isfield(synthetic_data.datasets(i).params, 'SNR')
            all_snr = [all_snr; synthetic_data.datasets(i).params.SNR];
        end
        % Handle both defect_density (new) and theta_cap (old) naming
        if isfield(synthetic_data.datasets(i).params, 'defect_density')
            all_theta = [all_theta; synthetic_data.datasets(i).params.defect_density];
        elseif isfield(synthetic_data.datasets(i).params, 'theta_cap')
            all_theta = [all_theta; synthetic_data.datasets(i).params.theta_cap];
        end
        % Handle both N_obs (new) and Nobs (old) naming
        if isfield(synthetic_data.datasets(i).params, 'N_obs')
            all_nobs = [all_nobs; synthetic_data.datasets(i).params.N_obs];
        elseif isfield(synthetic_data.datasets(i).params, 'Nobs')
            all_nobs = [all_nobs; synthetic_data.datasets(i).params.Nobs];
        end
        % Check for repetition parameter (check both 'rep' and 'repetition' for compatibility)
        if isfield(synthetic_data.datasets(i).params, 'rep')
            all_reps = [all_reps; synthetic_data.datasets(i).params.rep];
        elseif isfield(synthetic_data.datasets(i).params, 'repetition')
            all_reps = [all_reps; synthetic_data.datasets(i).params.repetition];
        end
    end
    
    % Get unique values
    unique_snr = unique(all_snr);
    unique_theta = unique(all_theta);
    unique_nobs = unique(all_nobs);
    
    % Determine repetition dimension size using two methods (robust approach)
    num_datasets = length(synthetic_data.datasets);
    if isfield(synthetic_data, 'param_sets') && ~isempty(synthetic_data.param_sets)
        num_unique_combinations = size(synthetic_data.param_sets, 1);
        % Method 2 (more robust): calculate from dataset count / unique combinations
        max_rep_calculated = num_datasets / num_unique_combinations;
        if mod(num_datasets, num_unique_combinations) == 0 && max_rep_calculated >= 1
            max_rep = round(max_rep_calculated);
        else
            % If not evenly divisible, use method 1 or default
            if ~isempty(all_reps)
                max_rep = max(all_reps);
            else
                max_rep = 1;
            end
        end
    else
        % Method 1: use rep field from datasets if available
        if ~isempty(all_reps)
            max_rep = max(all_reps);
        else
            max_rep = 1;  % No repetitions, use dimension 1 for backward compatibility
        end
    end
    
    unique_reps = 1:max_rep;  % Repetitions are 1, 2, 3, ..., max_rep
    
    % Store parameter values
    dataset_metrics.SNR_values = unique_snr;
    dataset_metrics.theta_cap_values = unique_theta;
    dataset_metrics.Nobs_values = unique_nobs;
    dataset_metrics.repetition_values = unique_reps;
    
    % Initialize metric arrays with pre-allocated memory (4D: SNR × theta × N_obs × rep)
    dims = [length(unique_snr), length(unique_theta), length(unique_nobs), max_rep];
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
            dataset_num = str2double(regexp(name, '\d+', 'match', 'once'));
            
            % Validate dataset number
            if isnan(dataset_num) || dataset_num < 1 || dataset_num > length(synthetic_data.datasets)
                warning('Invalid dataset number %d from file %s, skipping', dataset_num, result_files(i).name);
                continue;
            end
            
            % Get parameter values from dataset (handles repetitions)
            dataset = synthetic_data.datasets(dataset_num);
            
            % Extract parameters with backward compatibility
            if ~isfield(dataset.params, 'SNR')
                warning('Dataset %d missing SNR parameter, skipping', dataset_num);
                continue;
            end
            snr_val = dataset.params.SNR;
            
            % Handle both defect_density (new) and theta_cap (old) naming
            if isfield(dataset.params, 'defect_density')
                theta_val = dataset.params.defect_density;
            elseif isfield(dataset.params, 'theta_cap')
                theta_val = dataset.params.theta_cap;
            else
                warning('Dataset %d missing defect_density/theta_cap parameter, skipping', dataset_num);
                continue;
            end
            
            % Handle both N_obs (new) and Nobs (old) naming
            if isfield(dataset.params, 'N_obs')
                nobs_val = dataset.params.N_obs;
            elseif isfield(dataset.params, 'Nobs')
                nobs_val = dataset.params.Nobs;
            else
                warning('Dataset %d missing N_obs/Nobs parameter, skipping', dataset_num);
                continue;
            end
            
            % Get repetition index (check both 'rep' and 'repetition' fields, default to 1 for backward compatibility)
            if isfield(dataset.params, 'rep')
                rep_val = dataset.params.rep;
            elseif isfield(dataset.params, 'repetition')
                rep_val = dataset.params.repetition;
            else
                rep_val = 1;  % Backward compatibility: non-repetitive datasets use rep=1
            end
            
            % Map parameter values to indices (use exact matching for better accuracy)
            snr_idx = find(abs(unique_snr - snr_val) < 1e-10, 1);
            theta_idx = find(abs(unique_theta - theta_val) < 1e-10, 1);
            nobs_idx = find(abs(unique_nobs - nobs_val) < 1e-10, 1);
            rep_idx = rep_val;  % Repetition is already an index (1, 2, 3, ...)
            
            if isempty(snr_idx) || isempty(theta_idx) || isempty(nobs_idx)
                warning('Dataset %d: Could not map parameters to indices (SNR=%.2e, theta=%.2e, N_obs=%d), skipping', ...
                    dataset_num, snr_val, theta_val, nobs_val);
                continue;
            end
            
            % Validate repetition index
            if rep_idx < 1 || rep_idx > max_rep
                warning('Dataset %d: Invalid repetition index %d (max=%d), skipping', dataset_num, rep_idx, max_rep);
                continue;
            end
            
            indices = [snr_idx, theta_idx, nobs_idx, rep_idx];
            
            % Process current dataset
            try
                % Load result data
                data = load(fullfile(folder_path, result_files(i).name));
                
                % Store original and clean observations (4D indexing)
                dataset_metrics.Y{indices(1), indices(2), indices(3), indices(4)} = dataset.Y;
                dataset_metrics.Y_clean{indices(1), indices(2), indices(3), indices(4)} = dataset.Y_clean;
                dataset_metrics.A0_noiseless{indices(1), indices(2), indices(3), indices(4)} = dataset.A0_noiseless;
                % Store ground truth data
                if isfield(data, 'dataset_A0')
                    dataset_metrics.A0{indices(1), indices(2), indices(3), indices(4)} = data.dataset_A0{1};
                end
                if isfield(data, 'dataset_X0')
                    dataset_metrics.X0{indices(1), indices(2), indices(3), indices(4)} = data.dataset_X0{1};
                end
                
                % Store reconstruction data directly without flipping check
                if isfield(data, 'Xout')
                    dataset_metrics.Xout{indices(1), indices(2), indices(3), indices(4)} = data.Xout{1};
                end
                if isfield(data, 'Aout')
                    dataset_metrics.Aout{indices(1), indices(2), indices(3), indices(4)} = data.Aout{1};
                end
                if isfield(data, 'bout')
                    dataset_metrics.bout{indices(1), indices(2), indices(3), indices(4)} = data.bout{1};
                end
                
                % Calculate metrics for original order
                if isfield(data, 'Xout') && isfield(data, 'dataset_X0') && ...
                   ~isempty(data.Xout{1}) && ~isempty(data.dataset_X0{1}) && ...
                   isfield(data, 'Aout') && isfield(data, 'dataset_A0')
                    % Get quality metrics without flipping
                    %[~, quality_scores] = detect_kernel_flip(data.dataset_X0, data.Xout, data.dataset_A0, data.Aout, true);
                    kk = dataset.A0_noiseless;
                    [~, quality_scores] = detect_kernel_flip(data.dataset_X0, data.Xout, kk, data.Aout, true);
                    % Store metrics from original order
                    kernel_quality_final = mean(quality_scores.no_flip.kernel_similarity);
                    dataset_metrics.kernel_quality_final(indices(1), indices(2), indices(3), indices(4)) = kernel_quality_final;

                    act_sim_final = mean(quality_scores.no_flip.activation_similarity);
                    dataset_metrics.activation_similarity_final(indices(1), indices(2), indices(3), indices(4)) = act_sim_final;

                    % Calculate demixing score
                    [demixing_score, ~] = computeDemixingMetric(data.Xout{1});
                    dataset_metrics.demixing_score(indices(1), indices(2), indices(3), indices(4)) = demixing_score;
                    
                    % Calculate combined score
                    dataset_metrics.combined_activationScore(indices(1), indices(2), indices(3), indices(4)) = ...
                        computeCombined_activationScore(demixing_score, act_sim_final);
                end
                
                % Load metrics from trajectories
                extras = data.extras{1};
                dataset_metrics.extras{indices(1), indices(2), indices(3), indices(4)} = extras;
                if isstruct(extras) && isfield(extras, 'phase1')
                    % Store trajectories
                    if isfield(extras.phase1, 'kernel_quality_factors')
                        kq_traj = extras.phase1.kernel_quality_factors;
                        dataset_metrics.kernel_quality_trajectory{indices(1), indices(2), indices(3), indices(4)} = kq_traj;
                    end
                    
                    if isfield(extras.phase1, 'activation_metrics')
                        act_traj = extras.phase1.activation_metrics;
                        dataset_metrics.activation_similarity_trajectory{indices(1), indices(2), indices(3), indices(4)} = act_traj;
                    end
                    
                    % Store other phase1 metrics
                    if isfield(extras.phase1, 'residuals')
                        dataset_metrics.residuals{indices(1), indices(2), indices(3), indices(4)} = extras.phase1.residuals;
                    end
                    if isfield(extras.phase1, 'relative_changes')
                        dataset_metrics.relative_changes{indices(1), indices(2), indices(3), indices(4)} = extras.phase1.relative_changes;
                    end
                end
                
                % Store runtime
                if isfield(extras, 'runtime')
                    dataset_metrics.runtime(indices(1), indices(2), indices(3), indices(4)) = extras.runtime;
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
    
    % Show repetition dimension if it exists
    if isfield(metrics, 'repetition_values')
        fprintf('- Repetitions: %d values [%d to %d]\n', ...
            length(metrics.repetition_values), min(metrics.repetition_values), max(metrics.repetition_values));
        fprintf('- Array dimensions: %d × %d × %d × %d (SNR × theta × N_obs × rep)\n', ...
            length(metrics.SNR_values), length(metrics.theta_cap_values), ...
            length(metrics.Nobs_values), length(metrics.repetition_values));
    else
        fprintf('- Array dimensions: %d × %d × %d (SNR × theta × N_obs)\n', ...
            length(metrics.SNR_values), length(metrics.theta_cap_values), ...
            length(metrics.Nobs_values));
    end
    
    % Add summary of stored reconstruction data
    fprintf('\nReconstruction Data Summary:\n');
    non_empty_reconstructions = nnz(~cellfun(@isempty, metrics.Aout));
    total_points = numel(metrics.Aout);
    fprintf('- Stored reconstructions: %d/%d points (%.1f%%)\n', ...
        non_empty_reconstructions, total_points, ...
        100 * non_empty_reconstructions / total_points);
end 