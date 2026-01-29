function visualize_single_dataset()
    % Let user select the folder containing both parallel results and synthetic dataset
    default_path = fullfile(pwd, 'examples');
    folder_path = uigetdir(default_path, 'Select folder containing results');
    
    if folder_path == 0
        error('No folder selected');
    end
    
    % Find synthetic dataset file
    synthetic_files = dir(fullfile(folder_path, 'synthetic_datasets_*.mat'));
    if isempty(synthetic_files)
        error('No synthetic dataset file found in the selected folder');
    end
    
    % If multiple synthetic files exist, let user choose
    if length(synthetic_files) > 1
        [synthetic_file, ~] = uigetfile(fullfile(folder_path, 'synthetic_datasets_*.mat'), ...
            'Select synthetic dataset file');
        if isequal(synthetic_file, 0)
            error('No synthetic file selected');
        end
    else
        synthetic_file = synthetic_files(1).name;
    end
    
    % Load synthetic dataset
    synthetic_data = load(fullfile(folder_path, synthetic_file));
    
    % Get user input for dataset and parameter indices
    dataset_idx = input('Enter dataset index: ');
    param_idx = input('Enter parameter index: ');
    
    % Load parallel result file
    result_filename = sprintf('SBD_parallel_dataset%d_param_idx%d.mat', dataset_idx, param_idx);
    if ~exist(fullfile(folder_path, result_filename), 'file')
        error('Result file %s not found', result_filename);
    end
    
    result = load(fullfile(folder_path, result_filename));
    
    % Extract required data
    Y = synthetic_data.datasets(dataset_idx).Y;  % Original observation
    Y_clean = synthetic_data.datasets(dataset_idx).Y_clean;  % Clean observation
    A0 = result.dataset_A0{1};  % Ground truth kernels
    X0 = result.dataset_X0{1};  % Ground truth activations
    Aout = result.Aout{1};  % Reconstructed kernels
    Xout = result.Xout{1};  % Reconstructed activations
    bout = result.bout{1};  % Bias terms
    extras = result.extras{1};  % Convergence metrics
    
    % Calculate demixing metrics
    [demix_score, overlap_matrix] = computeDemixingMetric(Xout);
    extras.demixing_score = demix_score;
    extras.demixing_matrix = overlap_matrix;
    
    % Add clean observation to extras
    extras.Y_clean = Y_clean;
    
    % Call visualizeResults with dataset and parameter indices
    visualizeResults(Y, A0, Aout, X0, Xout, bout, extras, [dataset_idx, param_idx]);
    
    % Print dataset parameters
    param_combinations = result.param_combinations;
    lambda1 = param_combinations(param_idx, 1);
    mini_loop = param_combinations(param_idx, 2);
    
    fprintf('\nDataset %d Parameter Set %d:\n', dataset_idx, param_idx);
    fprintf('λ₁=%.3f, mini_loop=%d\n', lambda1, mini_loop);
    fprintf('theta_cap=%.0e, kernel_size=%.2f, SNR=%.1f\n', ...
        synthetic_data.datasets(dataset_idx).params.theta_cap, ...
        synthetic_data.datasets(dataset_idx).params.kernel_size(1), ...
        synthetic_data.datasets(dataset_idx).params.SNR);
    fprintf('Demixing Score: %.4f (higher is better)\n', demix_score);
end 