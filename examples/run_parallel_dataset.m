clc; clear;
run('../init_sbd');

%% Load pre-generated synthetic datasets
syn_loc = 'results\synthetic_datasets\synthetic_datasets_20251108_114222.mat';
load(syn_loc);  

%% Initialize kernels for all datasets
fprintf('Initializing kernels for %d datasets...\n', length(datasets));
% Choose initialization method ('activation' or 'random')
init_method = 'activation';
sigma_gaussian = 6.0;
if strcmp(init_method, 'activation')    
    A1_all = initialize_kernels_from_activations(datasets, sigma_gaussian);
else
    A1_all = initialize_kernels_random(datasets);
end

%% Save initialized kernels to synthetic dataset
for i = 1:length(datasets)
    datasets(i).A1_init = A1_all{i};
end

% Save updated datasets with initialized kernels
[filepath, name, ext] = fileparts(syn_loc);
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
maxIT = 15;

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
    parpool(10);  % Creates pool with 10 workers
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

% Create DataQueue for progress updates
D = parallel.pool.DataQueue;
completed_count = 0;
progress_callback = @(x) updateProgress(x, num_datasets);
afterEach(D, progress_callback);

% Parallel loop over datasets
parfor n = 1:num_datasets
    try
        % Extract dataset parameters
        dataset_Y = datasets(n).Y;
        dataset_kernel_size = datasets(n).params.kernel_size;
        dataset_A1 = A1_all{n};
        params_local = params_all{n};  % Get pre-processed params
        
        % Print start message
        fprintf('[START] Dataset %d/%d\n', n, num_datasets);

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

        % Save results
        filename = sprintf('SBD_parallel_dataset%d_optimal.mat', n);
        save_results(filename, Aout, Xout, bout, extras, param_combinations, ...
            param_idx, params_local.A0, params_local.X0);
        
        % Print end message and send progress update
        fprintf('[END] Dataset %d/%d (Runtime: %.2fs)\n', n, num_datasets, extras.runtime);
        send(D, n);
        
    catch ME
        fprintf('[ERROR] Dataset %d: %s\n', n, ME.message);
        send(D, n);  % Still count as completed for progress
    end
end

%% Cleanup worker configs
rmdir(config_dir, 's');

fprintf('Processing complete!\n');

%% Helper Functions
function A1_all = initialize_kernels_from_activations(datasets, sigma_gaussian)
    % Initialize kernels by finding most isolated activation points within allowed region
    % Falls back to closest-to-center if no valid points found
    % Input:
    %   datasets: 1×N struct array with fields Y, X0, params
    %   sigma_gaussian: Gaussian window parameter
    % Output:
    %   A1_all: 1×N cell array of initial kernel guesses
    
    num_datasets = length(datasets);
    A1_all = cell(1, num_datasets);
    isolation_threshold_factor = 10;  % Threshold = max/factor for defect detection
    
    for i = 1:num_datasets
        % Get data for this dataset
        Y = datasets(i).Y;
        X0 = datasets(i).X0;
        kernel_sizes = datasets(i).params.kernel_size;  % [K×2] matrix of kernel sizes
        num_kernels = size(X0, 3);
        
        % Get image dimensions
        img_height = size(Y, 1);
        img_width = size(Y, 2);
        
        % Initialize kernels for this dataset
        A1 = cell(1, num_kernels);
        most_isolated_points = cell(1, num_kernels);
        
        % Find defect positions for each kernel
        defect_positions = cell(1, num_kernels);
        for k = 1:num_kernels
            threshold = max(X0(:,:,k), [], 'all') / isolation_threshold_factor;
            [rows, cols] = find(X0(:,:,k) > threshold);
            defect_positions{k} = [rows, cols];
        end
        
        % Find most isolated points for each kernel
        for k = 1:num_kernels
            % Get kernel size
            ksize = kernel_sizes(k,:);
            half_size = floor(ksize / 2);
            
            % Define allowed region: centers where full kernel fits within image
            min_y = 1 + half_size(1);
            max_y = img_height - ksize(1) + 1 + half_size(1);
            min_x = 1 + half_size(2);
            max_x = img_width - ksize(2) + 1 + half_size(2);
            
            % Filter to points within allowed region
            valid_points = true(size(defect_positions{k}, 1), 1);
            for j = 1:size(defect_positions{k}, 1)
                y = defect_positions{k}(j,1);
                x = defect_positions{k}(j,2);
                if y < min_y || y > max_y || x < min_x || x > max_x
                    valid_points(j) = false;
                end
            end
            
            valid_defects = defect_positions{k}(valid_points,:);
            
            if isempty(valid_defects)
                % Fallback: use activation closest to image center
                if isempty(defect_positions{k})
                    warning('Dataset %d, Kernel %d: No activation points found. Using image center.', i, k);
                    img_center = [img_height/2, img_width/2];
                    most_isolated_points{k} = round(img_center);
                else
                    img_center = [img_height/2, img_width/2];
                    distances_to_center = sum((defect_positions{k} - img_center).^2, 2);
                    [~, closest_idx] = min(distances_to_center);
                    most_isolated_points{k} = defect_positions{k}(closest_idx,:);
                end
            else
                % Sum activation maps of all other kernels
                X_others = zeros(size(X0(:,:,1)));
                for l = 1:num_kernels
                    if l ~= k
                        X_others = X_others + X0(:,:,l);
                    end
                end
                
                % Get positions in other kernels
                other_threshold = max(X_others,[],'all') / isolation_threshold_factor;
                [other_rows, other_cols] = find(X_others > other_threshold);
                other_positions = [other_rows, other_cols];
                
                % Calculate isolation scores (min distance to other-kernel defects)
                S_k = zeros(size(valid_defects, 1), 1);
                for j = 1:size(valid_defects, 1)
                    diffs = other_positions - valid_defects(j,:);
                    distances = sum(diffs.^2, 2);
                    S_k(j) = min(distances);
                end
                
                % Select most isolated point
                [~, max_idx] = max(S_k);
                valid_indices = find(valid_points);
                max_idx = valid_indices(max_idx);
                most_isolated_points{k} = defect_positions{k}(max_idx,:);
            end
            
            % Extract kernel patch (pad with zeros if extends beyond boundaries)
            y_center = most_isolated_points{k}(1);
            x_center = most_isolated_points{k}(2);
            
            % Calculate extraction range (may extend beyond image)
            y_range_desired = (y_center - half_size(1)):(y_center + half_size(1));
            x_range_desired = (x_center - half_size(2)):(x_center + half_size(2));
            
            % Adjust for even-sized kernels
            if length(y_range_desired) > ksize(1)
                y_range_desired = y_range_desired(1:ksize(1));
            end
            if length(x_range_desired) > ksize(2)
                x_range_desired = x_range_desired(1:ksize(2));
            end
            
            % Clip to image boundaries
            y_range_valid = max(1, min(y_range_desired)):min(img_height, max(y_range_desired));
            x_range_valid = max(1, min(x_range_desired)):min(img_width, max(x_range_desired));
            
            % Extract and pad if needed
            if isempty(y_range_valid) || isempty(x_range_valid)
                warning('Dataset %d, Kernel %d: Entire kernel out of bounds. Using zeros.', i, k);
                kernel_patch = zeros(ksize(1), ksize(2));
            else
                kernel_patch_valid = Y(y_range_valid, x_range_valid);
                kernel_patch = zeros(ksize(1), ksize(2));
                
                % Place valid portion in correct position (accounting for clipping)
                y_offset = max(0, 1 - min(y_range_desired));
                x_offset = max(0, 1 - min(x_range_desired));
                y_kernel_start = y_offset + 1;
                y_kernel_end = y_offset + length(y_range_valid);
                x_kernel_start = x_offset + 1;
                x_kernel_end = x_offset + length(x_range_valid);
                
                kernel_patch(y_kernel_start:y_kernel_end, x_kernel_start:x_kernel_end) = kernel_patch_valid;
            end
            
            % Apply Gaussian window
            [h, w] = size(kernel_patch);
            [X_grid, Y_grid] = meshgrid(1:w, 1:h);
            center_x = (w + 1) / 2;
            center_y = (h + 1) / 2;
            gaussian_mask = exp(-((X_grid - center_x).^2 + (Y_grid - center_y).^2) / (2 * sigma_gaussian^2));
            A1{k} = kernel_patch .* gaussian_mask;
            
            % Normalize
            A1{k} = proj2oblique(A1{k});
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

%% Progress update function
function updateProgress(~, total_datasets)
    persistent count
    if isempty(count)
        count = 0;
    end
    count = count + 1;
    % Simple progress indicator - update every 5% or at completion
    update_interval = max(1, floor(total_datasets/20));
    if mod(count, update_interval) == 0 || count == total_datasets
        fprintf('Progress: %d/%d datasets completed (%.1f%%)\n', ...
            count, total_datasets, 100*count/total_datasets);
    end
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