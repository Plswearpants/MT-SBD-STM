%% Clear workspace and initialize
clc; clear;
run('../init_sbd');

%% 1. Initial Setup
% Load LDoS simulation data for kernel selection
load('example_data/LDoS_sim.mat');

% Display the 3D LDoS data for selection
fprintf('Displaying 3D LDoS simulation data...\n');
d3gridDisplay(LDoS_sim, 'dynamic');
title('LDoS Simulation Data - Use for Kernel Selection');

% Get user input for kernel slice selection
num_kernels = 2;  % Default number of kernels
fprintf('\nPlease select %d slice indices for kernels (1-%d):\n', num_kernels, size(LDoS_sim,3));

sliceidx = zeros(1, num_kernels);
for k = 1:num_kernels
    sliceidx(k) = input(sprintf('Enter slice index for kernel %d: ', k));
    while sliceidx(k) < 1 || sliceidx(k) > size(LDoS_sim,3)
        fprintf('Invalid slice index. Please enter a number between 1 and %d\n', size(LDoS_sim,3));
        sliceidx(k) = input(sprintf('Enter slice index for kernel %d: ', k));
    end
end

% Display and confirm kernel selection
confirm_kernel_selection(LDoS_sim, sliceidx);

% Define parameter ranges
theta_cap_values = logspace(-5, -3, 6);    % 6 different activation densities
area_ratio_values = linspace(0.01, 0.16, 6);  % 6 different kernel sizes
SNR_values = [10, 3.16, 1];                % 3 different noise levels

% Fixed parameters
fixed_params.num_kernels = num_kernels;
fixed_params.image_size = [500, 500];
fixed_params.relative_theta_cap = ones(1, num_kernels);
fixed_params.relative_kernel_size = ones(1, num_kernels);

%% 2. Generate Base Activations (Step 1)
fprintf('\nStep 1: Generating base activations for different theta_cap values...\n');
base_activations = struct('X0', {}, 'theta_cap', {});

% Create figure for activation review
act_fig = figure('Name', 'Base Activation Pattern Review', ...
                'Position', [100 100 800 800]);  % Made square for better visualization

% Generate and confirm base activations for each theta
for i = 1:length(theta_cap_values)
    confirmed = false;
    while ~confirmed
        fprintf('Generating activation for theta_cap = %.2e (%d/%d)\n', ...
            theta_cap_values(i), i, length(theta_cap_values));
        
        % Generate only activation maps (no kernels or noise yet)
        X0 = generate_activation_maps(theta_cap_values(i), fixed_params);
        
        % Display for review
        clf(act_fig);
        
        % Create RGB image combining both activations
        rgb_activation = zeros([fixed_params.image_size 3]);
        rgb_activation(:,:,1) = X0(:,:,1);  % First activation in red channel
        rgb_activation(:,:,2) = X0(:,:,2);  % Second activation in green channel
        
        % Display combined activations
        imagesc(rgb_activation);
        title(sprintf('Combined Activations (θ=%.2e)\nRed: Kernel 1, Green: Kernel 2\nYellow: Overlap', ...
            theta_cap_values(i)), 'FontSize', 12);
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
    base_activations(i).theta_cap = theta_cap_values(i);
end
close(act_fig);

%% 3. Create Kernel Size Variations (Step 2)
fprintf('\nStep 2: Creating kernel size variations...\n');
kernel_variations = cell(length(theta_cap_values), length(area_ratio_values));

for i = 1:length(theta_cap_values)
    for j = 1:length(area_ratio_values)
        fprintf('Processing area_ratio = %.3f for theta_cap = %.2e\n', ...
            area_ratio_values(j), theta_cap_values(i));
        
        % Use stored activation pattern
        X0 = base_activations(i).X0;
        
        % Generate kernels for this size
        [A0, A0_noiseless] = generate_kernels(LDoS_sim(:,:,sliceidx(1)), SNR_values(1), ...
            size(LDoS_sim(:,:,sliceidx(1)), 1), fixed_params.image_size(1), fixed_params.image_size(1)/size(LDoS_sim(:,:,sliceidx(1)), 1));
        
        % Create clean observation
        Y_clean = generate_clean_observation(A0_noiseless, X0);
        
        % Store without noise
        kernel_variations{i,j} = struct('X0', X0, ...
                                      'A0', {A0}, ...
                                      'A0_noiseless', {A0_noiseless}, ...
                                      'Y_clean', Y_clean);
    end
end

%% 4. Add Noise Variations (Step 3)
fprintf('\nGenerating noise variations...\n');

% Reverse SNR values to process low SNR (harder cases) first
SNR_values = flip(SNR_values);  % Now [1, 3.16, 10] instead of [10, 3.16, 1]

final_datasets = cell(length(theta_cap_values), length(area_ratio_values), length(SNR_values));
descriptions = cell(length(theta_cap_values), length(area_ratio_values), length(SNR_values));

% Reorder the loops to prioritize low SNR processing
for k = 1:length(SNR_values)
    for i = 1:length(theta_cap_values)
        for j = 1:length(area_ratio_values)
            fprintf('Processing SNR = %.1f, theta_cap = %.2e, area_ratio = %.3f (%d,%d,%d/%d,%d,%d)\n', ...
                SNR_values(k), theta_cap_values(i), area_ratio_values(j), ...
                k, i, j, length(SNR_values), length(theta_cap_values), length(area_ratio_values));
            
            % Get clean data from kernel variations
            clean_data = kernel_variations{i,j};
            
            % Add noise
            [Y, A0] = add_noise_to_dataset(clean_data.Y_clean, clean_data.A0_noiseless, SNR_values(k));
            
            % Store dataset
            final_datasets{i,j,k} = store_dataset(Y, clean_data.Y_clean, A0, ...
                clean_data.A0_noiseless, clean_data.X0, randn, ...
                theta_cap_values(i), area_ratio_values(j), SNR_values(k));
            
            % Create description
            descriptions{i,j,k} = sprintf('θ=%.2e, A_ratio=%.3f, SNR=%.1f', ...
                theta_cap_values(i), area_ratio_values(j), SNR_values(k));
        end
    end
end

% Create parameter set matrix (maintain the new SNR order)
[T, A, S] = meshgrid(theta_cap_values, area_ratio_values, SNR_values);
param_sets = [T(:), A(:), S(:)];  % Each row: [theta_cap, area_ratio, SNR]

%% 5. Save Results
% Convert cell arrays to 1D structure arrays with maintained order
datasets = convert_to_struct(final_datasets);  % Now returns 1D struct array
descriptions = reshape(descriptions, 1, []);  % Flatten to 1D

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
fprintf('Saved %d datasets as a 1D array\n', numel(datasets));
fprintf('Parameter space: %d theta_cap × %d area_ratio × %d SNR values\n', ...
    length(theta_cap_values), length(area_ratio_values), length(SNR_values));
fprintf('Parameter set matrix size: %d × 3\n', size(param_sets, 1));

%% 6. Display Results Overview
fprintf('\nDataset Overview:\n');
fprintf('Total number of datasets: %d\n\n', numel(datasets));

% Create 3D array of all observations
all_observations = zeros([fixed_params.image_size numel(datasets)]);

% Fill the 3D array in order: theta -> area_ratio -> SNR
for idx = 1:numel(datasets)
    all_observations(:,:,idx) = datasets(idx).Y;
end

% Display using d3gridDisplay
figure('Name', 'All Synthetic Observations');
d3gridDisplay(all_observations, 'dynamic');
title('All Synthetic Observations');

%% Helper Functions
function X0 = generate_activation_maps(N_obs, rho_d, p_scale, num_kernels)
    % Initialize output
    X0 = zeros(N_obs*p_scale, N_obs*p_scale, num_kernels);
    for k = 1: num_kernels
        % Generate activation map with strict density control
        target_defects = round(N_obs * N_obs * rho_d);  % Expected number of defects
        tolerance = 0.1;  % 10% tolerance
        min_defects = max(round(target_defects * (1-tolerance)), 1);
        max_defects = round(target_defects * (1+tolerance));
        
        X_good = false;
        while ~X_good
            X = double(rand(N_obs, N_obs) <= rho_d);
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
    
    % Add noise to kernels
    avg_var = 0;
    for k = 1:length(A0_noiseless)
        avg_var = avg_var + var(A0_noiseless{k}, 0, "all");
    end
    avg_var = avg_var / length(A0_noiseless);
    
    eta_kernel = avg_var/SNR;
    for k = 1:length(A0_noiseless)
        A0{k} = A0_noiseless{k} + sqrt(eta_kernel)*randn(size(A0_noiseless{k}));
        A0{k} = proj2oblique(A0{k});
    end
    
    % Add noise to observation
    eta = var(Y_clean, 0, "all") / SNR;
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

function dataset = store_dataset(Y, Y_clean, A0, A0_noiseless, X0, b0, theta_cap, area_ratio, SNR)
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
    
    dataset.params = struct('theta_cap', theta_cap, ...
                          'kernel_size', kernel_sizes, ...  % [2×2] matrix of kernel sizes
                          'SNR', SNR);
    dataset.b0 = b0;
end

function datasets_struct = convert_to_struct(datasets_cell)
    % Convert cell array to 1D struct array
    datasets_struct = [datasets_cell{:}];  % Flatten to 1D
    datasets_struct = reshape(datasets_struct, 1, []);  % Ensure its a row vector
end

function regenerate_dataset()
    uiresume;
end

function confirm_dataset()
    assignin('base', 'confirmed', true);
    uiresume;
end 