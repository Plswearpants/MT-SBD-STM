function dataset_metrics = load_datasets_metrics()
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
        param_sets = synthetic_data.param_sets;  % [theta, Area_ratio, SNR]
    catch ME
        error('Failed to load synthetic dataset file: %s', ME.message);
    end

    % Create parameter mapping (hash table approach)
    [unique_theta, ~, theta_idx] = unique(param_sets(:,1));
    [unique_area, ~, area_idx] = unique(param_sets(:,2));
    [unique_snr, ~, snr_idx] = unique(param_sets(:,3));
    
    % Store parameter values
    dataset_metrics.theta_cap_values = unique_theta;
    dataset_metrics.area_ratio_values = unique_area;
    dataset_metrics.SNR_values = unique_snr;
    
    % Create parameter lookup table
    param_lookup = zeros(size(param_sets, 1), 3);
    for i = 1:size(param_sets, 1)
        param_lookup(i,:) = [theta_idx(i), area_idx(i), snr_idx(i)];
    end
    
    % Initialize metric arrays with pre-allocated memory
    dims = [length(unique_theta), length(unique_area), length(unique_snr)];
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
                
                % Store ground truth data
                if isfield(data, 'dataset_A0')
                    dataset_metrics.A0{indices(1), indices(2), indices(3)} = data.dataset_A0{1};
                end
                if isfield(data, 'dataset_X0')
                    dataset_metrics.X0{indices(1), indices(2), indices(3)} = data.dataset_X0{1};
                end
                
                % Handle kernel flipping if necessary
                is_flipped = false;
                if isfield(data, 'Xout') && isfield(data, 'dataset_X0') && ...
                   ~isempty(data.Xout{1}) && ~isempty(data.dataset_X0{1}) && ...
                   isfield(data, 'Aout') && isfield(data, 'dataset_A0')
                    % Check if the kernels are flipped
                    [is_flipped, quality_scores, dataset_metrics.Xout{indices(1), indices(2), indices(3)}, dataset_metrics.Aout{indices(1), indices(2), indices(3)}] = detect_kernel_flip(data.dataset_X0, data.Xout, data.dataset_A0, data.Aout);
                    if is_flipped
                        fprintf('Dataset %d: Kernels flipped, correcting order...\n', dataset_num);

                        % recalculate metrics for flipped kernels and activation
                        
                        kernel_quality_final = quality_scores.flipped.kernel_similarity;
                        dataset_metrics.kernel_quality_final(indices(1), indices(2), indices(3)) = kernel_quality_final;

                        act_sim_final = quality_scores.flipped.activation_similarity;
                        dataset_metrics.activation_similarity_final(indices(1), indices(2), indices(3)) = act_sim_final;

                        % demixing score
                        [demixing_score, ~] = computeDemixingMetric(dataset_metrics.Xout{indices(1), indices(2), indices(3)});
                        dataset_metrics.demixing_score(indices(1), indices(2), indices(3)) = demixing_score;
                        
                        % combined score
                        dataset_metrics.combined_activationScore(indices(1), indices(2), indices(3)) = ...
                            computeCombined_activationScore(demixing_score, act_sim_final);
                    end
                end
                
                % Store reconstruction data
                if isfield(data, 'Aout')
                    dataset_metrics.Aout{indices(1), indices(2), indices(3)} = data.Aout{1};
                end
                if isfield(data, 'Xout')
                    dataset_metrics.Xout{indices(1), indices(2), indices(3)} = data.Xout{1};
                end
                if isfield(data, 'bout')
                    dataset_metrics.bout{indices(1), indices(2), indices(3)} = data.bout{1};
                end
                
                % Calculate or load metrics
                if ~is_flipped
                    % Load metrics from trajectories
                    extras = data.extras{1};
                    dataset_metrics.extras{indices(1), indices(2), indices(3)} = extras;
                    if isstruct(extras) && isfield(extras, 'phase1')
                        % Store trajectories
                        if isfield(extras.phase1, 'kernel_quality_factors')
                            kq_traj = extras.phase1.kernel_quality_factors;
                            dataset_metrics.kernel_quality_trajectory{indices(1), indices(2), indices(3)} = kq_traj;
                            if ~isempty(kq_traj)
                                dataset_metrics.kernel_quality_final(indices(1), indices(2), indices(3)) = prod(kq_traj(end,:));
                            end
                        end
                        
                        if isfield(extras.phase1, 'activation_metrics')
                            act_traj = extras.phase1.activation_metrics;
                            dataset_metrics.activation_similarity_trajectory{indices(1), indices(2), indices(3)} = act_traj;
                            if ~isempty(act_traj)
                                dataset_metrics.activation_similarity_final(indices(1), indices(2), indices(3)) = act_traj(end);
                            end
                        end
                        
                        % Store other phase1 metrics
                        if isfield(extras.phase1, 'residuals')
                            dataset_metrics.residuals{indices(1), indices(2), indices(3)} = extras.phase1.residuals;
                        end
                        if isfield(extras.phase1, 'relative_changes')
                            dataset_metrics.relative_changes{indices(1), indices(2), indices(3)} = extras.phase1.relative_changes;
                        end
                        
                        % Compute demixing score from final activation maps
                        if isfield(data, 'Xout') && ~isempty(data.Xout{1})
                            [demixing_score, ~] = computeDemixingMetric(data.Xout{1});
                            dataset_metrics.demixing_score(indices(1), indices(2), indices(3)) = demixing_score;
                        end

                        % Calculate combined score for non-flipped case
                        d_final = dataset_metrics.demixing_score(indices(1), indices(2), indices(3));
                        a_final = dataset_metrics.activation_similarity_final(indices(1), indices(2), indices(3));
                        if ~isnan(d_final) && ~isnan(a_final)
                            dataset_metrics.combined_activationScore(indices(1), indices(2), indices(3)) = ...
                                computeCombined_activationScore(d_final, a_final);
                        end
                    end
                    
                    % Store runtime
                    if isfield(extras, 'runtime')
                        dataset_metrics.runtime(indices(1), indices(2), indices(3)) = extras.runtime;
                    end
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
    fprintf('- Theta cap: %d values [%.2e to %.2e]\n', ...
        length(metrics.theta_cap_values), min(metrics.theta_cap_values), max(metrics.theta_cap_values));
    fprintf('- Area ratio: %d values [%.2f to %.2f]\n', ...
        length(metrics.area_ratio_values), min(metrics.area_ratio_values), max(metrics.area_ratio_values));
    fprintf('- SNR: %d values [%.2e to %.2e]\n', ...
        length(metrics.SNR_values), min(metrics.SNR_values), max(metrics.SNR_values));
    
    % Add summary of stored reconstruction data
    fprintf('\nReconstruction Data Summary:\n');
    non_empty_reconstructions = nnz(~cellfun(@isempty, metrics.Aout));
    total_points = numel(metrics.Aout);
    fprintf('- Stored reconstructions: %d/%d points (%.1f%%)\n', ...
        non_empty_reconstructions, total_points, ...
        100 * non_empty_reconstructions / total_points);
end 