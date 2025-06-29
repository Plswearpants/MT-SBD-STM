%% Clear workspace and initialize
clc; clear;
run('../init_sbd');

%% 1. Initial Setup
% Load LDoS simulation data for kernel selection
load('example_data/LDoS_single_defect_save.mat');
LDoS_sim = LDoS_result;
% Display the 3D LDoS data for selection
fprintf('Displaying 3D LDoS simulation data...\n');
figure;
d3gridDisplay(LDoS_sim, 'dynamic');
title('LDoS Simulation Data - Use for Kernel Selection');

% Get user input for kernel slice selection
fprintf('\nPlease select slice indices (space-separated numbers): ');
input_str = input('', 's');
sliceidx = str2num(input_str);
num_kernels = length(sliceidx);

% Fixed parameters
fixed_params.p_scale = 3;                       % Resolution factor
fixed_params.N_single = N;                      % Input lattice size
fixed_params.num_kernels = num_kernels;         % Number of kernels

% Define parameter ranges for param_sets
SNR_values = [3, 2, 1];                     % Different noise levels
defect_density_values = logspace(-3, -1.5, 5);    % Different activation densities
N_obs_values = [50, 100, 150, 200];           % Different observation lattice sizes

% Create parameter set matrix
[S, D, N] = meshgrid(SNR_values, defect_density_values, N_obs_values);
param_sets = [S(:), D(:), N(:)];  % Each row: [SNR, defect_density, N_obs]

fprintf('Parameter space setup complete:\n');
fprintf('- SNR values: %d points from %.2f to %.2f\n', length(SNR_values), min(SNR_values), max(SNR_values));
fprintf('- Defect density values: %d points from %.2e to %.2e\n', length(defect_density_values), min(defect_density_values), max(defect_density_values));
fprintf('- N_obs values: %d points from %d to %d\n', length(N_obs_values), min(N_obs_values), max(N_obs_values));
fprintf('Total parameter combinations: %d\n', size(param_sets, 1));

% Display and confirm kernel selection
confirm_kernel_selection(LDoS_sim, sliceidx);

%% 2. Generate Base Activations (Step 1)
fprintf('\nStep 1: Generating base activations for different defect densities...\n');
base_activations = struct('X0', {}, 'defect_density', {}, 'N_obs', {});

% Create figure for activation review
act_fig = figure('Name', 'Base Activation Pattern Review', ...
                'Position', [100 100 800 800]);

% Get unique combinations of defect_density and N_obs
unique_combinations = unique(param_sets(:,[2 3]), 'rows');  % [rho_d, N_obs]

% Generate and confirm base activations for each combination
for i = 1:size(unique_combinations, 1)
    confirmed = false;
    while ~confirmed
        rho_d = unique_combinations(i,1);
        N_obs = unique_combinations(i,2);
        
        fprintf('Generating activation for rho_d = %.2e, N_obs = %d (%d/%d)\n', ...
            rho_d, N_obs, i, size(unique_combinations,1));
        
        % Generate activation maps for this N_obs and density
        X0 = generate_activation_maps(N_obs, rho_d, fixed_params.p_scale, fixed_params.num_kernels);
        
        % Display for review
        clf(act_fig);
        
        % Create RGB image combining both activations
        rgb_activation = zeros([N_obs*fixed_params.p_scale N_obs*fixed_params.p_scale 3]);
        rgb_activation(:,:,1) = X0(:,:,1);  % First activation in red channel
        rgb_activation(:,:,2) = X0(:,:,2);  % Second activation in green channel
        
        % Display combined activations
        imagesc(rgb_activation);
        title(sprintf('Combined Activations (ρ_d=%.2e, N_{obs}=%d)\nRed: Kernel 1, Green: Kernel 2\nYellow: Overlap', ...
            rho_d, N_obs), 'FontSize', 12);
        axis square;
        
        % Add colorbar with custom labels
        colorbar('Ticks', [0, 0.5, 1], ...
                'TickLabels', {'None', '', 'Active'});
        
        % Add UI controls
        regenerate_btn = uicontrol('Parent', act_fig, ...
                                 'Style', 'pushbutton', ...
                                 'String', 'Regenerate', ...
                                 'Position', [20 20 100 40], ...
                                 'Callback', @(src,event) regenerate_dataset());
        
        confirm_btn = uicontrol('Parent', act_fig, ...
                               'Style', 'pushbutton', ...
                               'String', 'Confirm', ...
                               'Position', [140 20 100 40], ...
                               'Callback', @(src,event) confirm_dataset());
        
        % Add text showing activation statistics
        stats_str = sprintf('Activation Rates:\nKernel 1: %.4f%%\nKernel 2: %.4f%%\nOverlap: %.4f%%', ...
            100*mean(X0(:,:,1),'all'), ...
            100*mean(X0(:,:,2),'all'), ...
            100*mean(X0(:,:,1) & X0(:,:,2), 'all'));
        
        uicontrol('Parent', act_fig, ...
                 'Style', 'text', ...
                 'String', stats_str, ...
                 'Position', [260 20 200 60], ...
                 'BackgroundColor', get(act_fig, 'Color'), ...
                 'HorizontalAlignment', 'left');
        
        drawnow;
        uiwait(act_fig);
    end
    
    % Store confirmed activation
    base_activations(i).X0 = X0;
    base_activations(i).defect_density = rho_d;
    base_activations(i).N_obs = N_obs;
end
close(act_fig);

%% 3. Create Dataset Variations
fprintf('\nGenerating datasets for all parameter combinations...\n');
final_datasets = cell(size(param_sets, 1), 1);
descriptions = cell(size(param_sets, 1), 1);

for i = 1:size(param_sets, 1)
    % Get parameters for this iteration
    SNR = param_sets(i,1);
    rho_d = param_sets(i,2);
    N_obs = param_sets(i,3);
    
    % Calculate M for this SNR
    [M, ~] = find_cutoff_noise_intersection(LDoS_sim(:,:,sliceidx(1)), SNR, fixed_params.N_single, false);
    
    % Calculate area_ratio after we have M and N_obs
    area_ratio = M / N_obs;
    
    fprintf('Processing combination %d/%d: SNR=%.1f, rho_d=%.2e, N_obs=%d (area_ratio=%.2e)\n', ...
        i, size(param_sets,1), SNR, rho_d, N_obs, area_ratio);
    
    % Find corresponding base activation
    base_idx = find([base_activations.defect_density] == rho_d & [base_activations.N_obs] == N_obs);
    X0 = base_activations(base_idx).X0;
    
    % Generate kernels
    [A0, A0_noiseless] = generate_kernels(LDoS_sim(:,:,sliceidx), SNR, ...
        fixed_params.N_single, N_obs, fixed_params.p_scale);
    
    % Generate clean observation
    Y_clean = generate_clean_observation(A0_noiseless, X0);
    
    % Add noise
    [Y, A0] = add_noise_to_dataset(Y_clean, A0_noiseless, SNR);
    
    % Store dataset with both N_obs and derived area_ratio
    final_datasets{i} = store_dataset(Y, Y_clean, A0, A0_noiseless, X0, randn, rho_d, SNR, N_obs, area_ratio);
    descriptions{i} = sprintf('SNR=%.1f, ρ_d=%.2e, N_obs=%d (area_ratio=%.2e)', ...
        SNR, rho_d, N_obs, area_ratio);
end

% Convert cell arrays to 1D structure arrays
datasets = [final_datasets{:}];
descriptions = reshape(descriptions, 1, []);

%% Save Results
% Add a note about the ordering in the saved file
ordering_info = struct('SNR_order', 'ascending', ...  % lowest to highest SNR
                      'SNR_values', SNR_values);

% Save results with timestamp and ordering information
timestamp = datestr(now, 'yyyymmdd_HHMMSS');
save_dir = 'results/synthetic_datasets';
if ~exist(save_dir, 'dir')
    mkdir(save_dir);
end

save_filename = fullfile(save_dir, sprintf('synthetic_datasets_%s.mat', timestamp));
save(save_filename, 'datasets', 'descriptions', 'param_sets', ...
    'fixed_params', 'sliceidx', 'LDoS_sim', 'ordering_info');

fprintf('\nDatasets saved to: %s\n', save_filename);
fprintf('Saved %d datasets\n', numel(datasets));
fprintf('Parameter combinations: %d (SNR × defect_density × N_obs)\n', size(param_sets, 1));

%% 6. Display Results Overview
fprintf('\nDataset Overview:\n');
fprintf('Total number of datasets: %d\n\n', numel(datasets));

% Get unique values for each parameter using arrayfun
unique_SNR = unique(arrayfun(@(x) x.params.SNR, datasets));
unique_rho_d = unique(arrayfun(@(x) x.params.defect_density, datasets));
unique_N_obs = unique(arrayfun(@(x) x.params.N_obs, datasets));

fprintf('Parameter ranges:\n');
fprintf('- SNR: %d values (%.1f to %.1f)\n', length(unique_SNR), min(unique_SNR), max(unique_SNR));
fprintf('- Defect density: %d values (%.2e to %.2e)\n', length(unique_rho_d), min(unique_rho_d), max(unique_rho_d));
fprintf('- N_obs: %d values (%d to %d)\n', length(unique_N_obs), min(unique_N_obs), max(unique_N_obs));

% Create separate figures for different N_obs values
for n = 1:length(unique_N_obs)
    % Get datasets for this N_obs
    idx_N_obs = find(arrayfun(@(x) x.params.N_obs == unique_N_obs(n), datasets));
    N_current = length(idx_N_obs);
    
    % Create 3D array for this N_obs
    obs_size = size(datasets(idx_N_obs(1)).Y);
    all_obs_N = zeros([obs_size N_current]);
    
    for i = 1:N_current
        all_obs_N(:,:,i) = datasets(idx_N_obs(i)).Y;
    end
    
    % Create new figure for this N_obs
    figure('Name', sprintf('Synthetic Observations N_obs=%d', unique_N_obs(n)));
    d3gridDisplay(all_obs_N, 'dynamic');
    title(sprintf('N_{obs} = %d (%d datasets)', unique_N_obs(n), N_current));
end

% Add a summary of area ratios
figure('Name', 'Area Ratio Distribution');
area_ratios = arrayfun(@(x) x.params.area_ratio, datasets);
histogram(area_ratios, 20);
title('Distribution of Area Ratios');
xlabel('Area Ratio (M/N_{obs})');
ylabel('Count');

% Print summary statistics
fprintf('\nArea Ratio Statistics:\n');
fprintf('- Min: %.2e\n', min(area_ratios));
fprintf('- Max: %.2e\n', max(area_ratios));
fprintf('- Mean: %.2e\n', mean(area_ratios));
fprintf('- Median: %.2e\n', median(area_ratios));

%% Helper Functions
function X0 = generate_activation_maps(N_single, rho_d, p_scale, num_kernels)
    % Initialize output
    X0 = zeros(N_single*p_scale, N_single*p_scale, num_kernels);
    for k = 1: num_kernels
        % Generate activation map with strict density control
        target_defects = round(N_single * N_single * rho_d);  % Expected number of defects
        tolerance = 0.1;  % 10% tolerance
        min_defects = max(round(target_defects * (1-tolerance)), 1);
        max_defects = round(target_defects * (1+tolerance));
        
        X_good = false;
        while ~X_good
            X = double(rand(N_single, N_single) <= rho_d);
            num_defects = sum(X(:));
            X_good = (num_defects >= min_defects) && (num_defects <= max_defects);
        end
        
        % Upsample X
        X0(:,:,k) = upsample_with_zero_blocks(X, p_scale);
    end
end

function [A0, A0_noiseless] = generate_kernels(rho_single, SNR, N_single, N_obs, p_scale)
    % Generate kernels using SNR-based cutoff and proper scaling
    % 
    % Inputs:
    %   rho_single: Single defect QPI simulation (3D array for multiple slices)
    %   SNR: Signal-to-noise ratio
    %   N_single: Lattice size of input QPI pattern
    %   N_obs: Size of observation lattice
    %   p_scale: Resolution factor (pixels per lattice site)
    %
    % Outputs:
    %   A0: Cell array of noisy kernels
    %   A0_noiseless: Cell array of noiseless kernels
    
    num_kernels = size(rho_single, 3);
    A0 = cell(1, num_kernels);
    A0_noiseless = cell(1, num_kernels);
    
    % Process each kernel
    for k = 1:num_kernels
        % Find cutoff M for this slice
        [M, M_pixels, ~] = find_cutoff_noise_intersection(rho_single(:,:,k), SNR, N_single, false);
        
        % Get center and crop rho_single
        [ny, nx] = size(rho_single(:,:,k));
        center = ceil([ny, nx]/2);
        
        % Crop the pattern
        range_y = max(1, center(1)-M_pixels):min(ny, center(1)+M_pixels);
        range_x = max(1, center(2)-M_pixels):min(nx, center(2)+M_pixels);
        rho_single_cutoff = rho_single(range_y, range_x, k);
        
        % Resize cutoff pattern to match observation scale
        target_size = round([2*M*p_scale 2*M*p_scale]);
        A0_noiseless{k} = imresize(rho_single_cutoff, target_size);
        A0_noiseless{k} = proj2oblique(A0_noiseless{k});
        
        % Add noise based on SNR
        signal_variance = var(A0_noiseless{k}(:));
        eta = signal_variance / SNR;
        A0{k} = A0_noiseless{k} + sqrt(eta) * randn(size(A0_noiseless{k}));
        A0{k} = proj2oblique(A0{k});
    end
end

function [Y, A0] = add_noise_to_dataset(Y_clean, A0_noiseless, SNR)
    % Add noise to clean observation and kernels
    % Initialize noisy kernels
    A0 = cell(size(A0_noiseless));
    
    % Calculate mean variance of kernels for noise level determination
    kernel_variances = zeros(1, length(A0_noiseless));
    for k = 1:length(A0_noiseless)
        kernel_variances(k) = var(A0_noiseless{k}, [], 'all');
    end
    mean_kernel_variance = mean(kernel_variances);
    
    % Calculate noise variance based on kernel variance
    eta = mean_kernel_variance / SNR;
    
    % Add noise to kernels
    for k = 1:length(A0_noiseless)
        A0{k} = A0_noiseless{k} + sqrt(eta) * randn(size(A0_noiseless{k}));
        A0{k} = proj2oblique(A0{k});
    end
    
    % Add noise to observation using the same eta based on kernel variance
    Y = Y_clean + sqrt(eta) * randn(size(Y_clean));
end

function Y_clean = generate_clean_observation(A0_noiseless, X0)
    % Generate clean observation by convolving kernels with activation maps
    Y_clean = zeros(size(X0, 1:2));
    
    for k = 1:size(X0, 3)
        Y_clean = Y_clean + convfft2(A0_noiseless{k}, X0(:,:,k));
    end
end

function confirm_kernel_selection(LDoS_sim, sliceidx)
    figure('Name', 'Selected Kernel Slices');
    for k = 1:length(sliceidx)
        subplot(1, length(sliceidx), k);
        imagesc(LDoS_sim(:,:,sliceidx(k)));
        title(sprintf('Kernel %d (Slice %d)', k, sliceidx(k)));
        colorbar;
        axis square;
    end
    
    confirmation = input('\nAre you satisfied with these kernel selections? (y/n): ', 's');
    while ~strcmpi(confirmation, 'y')
        % Allow reselection if not satisfied
        fprintf('\nPlease reselect slices:\n');
        for k = 1:length(sliceidx)
            sliceidx(k) = input(sprintf('Enter slice index for kernel %d: ', k));
            while sliceidx(k) < 1 || sliceidx(k) > size(LDoS_sim,3)
                fprintf('Invalid slice index. Please enter a number between 1 and %d\n', size(LDoS_sim,3));
                sliceidx(k) = input(sprintf('Enter slice index for kernel %d: ', k));
            end
        end
        
        % Update display
        for k = 1:length(sliceidx)
            subplot(1, length(sliceidx), k);
            imagesc(LDoS_sim(:,:,sliceidx(k)));
            title(sprintf('Kernel %d (Slice %d)', k, sliceidx(k)));
            colorbar;
            axis square;
        end
        
        confirmation = input('\nAre you satisfied with these kernel selections? (y/n): ', 's');
    end
end

function dataset = store_dataset(Y, Y_clean, A0, A0_noiseless, X0, b0, rho_d, SNR, N_obs, area_ratio)
    dataset.Y = Y;
    dataset.Y_clean = Y_clean;
    dataset.A0 = A0;
    dataset.A0_noiseless = A0_noiseless;
    dataset.X0 = X0;
    
    % Calculate kernel sizes as 2×2 matrix where each row is [height, width] for one kernel
    kernel_sizes = zeros(2, 2);  % Initialize [2×2] matrix
    for k = 1:length(A0)
        kernel_sizes(k,:) = size(A0{k});  % Store [height, width] for each kernel
    end
    
    % Initialize A1 using most isolated points method
    dataset.A1 = initialize_kernels_from_isolated_points(Y, X0, A0);
    
    dataset.params = struct('defect_density', rho_d, ...
                          'kernel_size', kernel_sizes, ...  % [n×2] matrix of kernel sizes
                          'SNR', SNR, ...
                          'N_obs', N_obs, ...
                          'area_ratio', area_ratio);
    dataset.b0 = b0;
end

function A1 = initialize_kernels_from_isolated_points(Y, X0, A0)
    % Initialize kernels by finding most isolated activation points
    % Inputs:
    %   Y: Observation data
    %   X0: Ground truth activation maps (binary)
    %   A0: Ground truth kernels (for size reference)
    % Output:
    %   A1: Cell array of initialized kernels with same sizes as A0
    
    num_kernels = size(X0, 3);
    A1 = cell(1, num_kernels);
    
    % Get kernel sizes from A0
    kernel_sizes = zeros(num_kernels, 2);
    for k = 1:num_kernels
        kernel_sizes(k,:) = size(A0{k});
    end
    
    % Find most isolated points for each kernel
    most_isolated_points = cell(1, num_kernels);
    
    for k = 1:num_kernels
        % Get positions of defects (all activations are 1 in synthetic data)
        [rows, cols] = find(X0(:,:,k));
        defect_positions = [rows, cols];
        
        if isempty(defect_positions)
            warning('No defects found for kernel %d', k);
            A1{k} = zeros(kernel_sizes(k,:));  % Initialize with zeros of correct size
            continue;
        end
        
        % Create summed activation map of all other kernels
        X_others = zeros(size(X0(:,:,1)));
        for l = 1:num_kernels
            if l ~= k
                X_others = X_others + X0(:,:,l);
            end
        end
        
        % Get positions of defects in other kernels
        [other_rows, other_cols] = find(X_others);
        other_positions = [other_rows, other_cols];
        
        % Calculate isolation scores for all points
        S_k = zeros(size(defect_positions, 1), 1);
        for i = 1:size(defect_positions, 1)
            diffs = other_positions - defect_positions(i,:);
            distances = sum(diffs.^2, 2);
            S_k(i) = min(distances);
        end
        
        % Get exact kernel size from A0
        [kernel_h, kernel_w] = size(A0{k});
        half_h = floor(kernel_h/2);
        half_w = floor(kernel_w/2);
        
        % Filter points that are away from boundaries
        valid_points = true(size(defect_positions, 1), 1);
        for i = 1:size(defect_positions, 1)
            y = defect_positions(i,1);
            x = defect_positions(i,2);
            
            if y <= half_h || y >= size(X0,1) - half_h || ...
               x <= half_w || x >= size(X0,2) - half_w
                valid_points(i) = false;
            end
        end
        
        % If we have valid points away from boundaries, use the one with highest score
        if any(valid_points)
            valid_defects = defect_positions(valid_points,:);
            valid_scores = S_k(valid_points);
            [~, max_idx] = max(valid_scores);
            most_isolated_points{k} = valid_defects(max_idx,:);
            
            % Extract kernel directly since we know it's away from boundaries
            y = most_isolated_points{k}(1);
            x = most_isolated_points{k}(2);
            A1{k} = Y(y-half_h:y+half_h, x-half_w:x+half_w);
            
            % Ensure exact size match with A0
            if ~isequal(size(A1{k}), size(A0{k}))
                A1{k} = imresize(A1{k}, size(A0{k}));
            end
        else
            % Fallback: use the point with highest isolation score and pad
            [~, max_idx] = max(S_k);
            most_isolated_points{k} = defect_positions(max_idx,:);
            
            % Extract kernel with padding
            y = most_isolated_points{k}(1);
            x = most_isolated_points{k}(2);
            
            % Calculate valid ranges that stay within image boundaries
            y_start = max(1, y - half_h);
            y_end = min(size(Y,1), y + half_h);
            x_start = max(1, x - half_w);
            x_end = min(size(Y,2), x + half_w);
            
            % Extract the valid portion
            kernel_patch = Y(y_start:y_end, x_start:x_end);
            
            % Create zero-padded kernel of required size
            A1{k} = zeros(kernel_sizes(k,:));
            
            % Calculate offsets for centering the patch
            y_offset = half_h - (y - y_start);
            x_offset = half_w - (x - x_start);
            
            % Place the patch in the center of the zero-padded kernel
            A1{k}(y_offset + (1:size(kernel_patch,1)), x_offset + (1:size(kernel_patch,2))) = kernel_patch;
            
            % Ensure exact size match with A0
            if ~isequal(size(A1{k}), size(A0{k}))
                A1{k} = imresize(A1{k}, size(A0{k}));
            end
        end
        
        % Normalize the kernel
        A1{k} = proj2oblique(A1{k});
        
        % Final size check
        assert(isequal(size(A1{k}), size(A0{k})), ...
            'Kernel size mismatch: A1{%d} size %dx%d does not match A0{%d} size %dx%d', ...
            k, size(A1{k},1), size(A1{k},2), k, size(A0{k},1), size(A0{k},2));
    end
end

function regenerate_dataset()
    uiresume;
end

function confirm_dataset()
    assignin('base', 'confirmed', true);
    uiresume;
end 