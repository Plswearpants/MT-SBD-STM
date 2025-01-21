%% Clear workspace and initialize
clc; clear;
run('../init_sbd');  % Adjust path as needed

%% Load and display LDoS simulation data for kernel selection
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
    
    % Validate input
    while sliceidx(k) < 1 || sliceidx(k) > size(LDoS_sim,3)
        fprintf('Invalid slice index. Please enter a number between 1 and %d\n', size(LDoS_sim,3));
        sliceidx(k) = input(sprintf('Enter slice index for kernel %d: ', k));
    end
end

% Display selected slices for confirmation
figure('Name', 'Selected Kernel Slices');
for k = 1:num_kernels
    subplot(1, num_kernels, k);
    imagesc(LDoS_sim(:,:,sliceidx(k)));
    title(sprintf('Kernel %d (Slice %d)', k, sliceidx(k)));
    colorbar;
    axis square;
end

% Confirm selection
confirmation = input('\nAre you satisfied with these kernel selections? (y/n): ', 's');
while ~strcmpi(confirmation, 'y')
    % If not satisfied, allow reselection
    fprintf('\nPlease reselect slices:\n');
    for k = 1:num_kernels
        sliceidx(k) = input(sprintf('Enter slice index for kernel %d: ', k));
        
        % Validate input
        while sliceidx(k) < 1 || sliceidx(k) > size(LDoS_sim,3)
            fprintf('Invalid slice index. Please enter a number between 1 and %d\n', size(LDoS_sim,3));
            sliceidx(k) = input(sprintf('Enter slice index for kernel %d: ', k));
        end
    end
    
    % Update display
    for k = 1:num_kernels
        subplot(1, num_kernels, k);
        imagesc(LDoS_sim(:,:,sliceidx(k)));
        title(sprintf('Kernel %d (Slice %d)', k, sliceidx(k)));
        colorbar;
        axis square;
    end
    
    confirmation = input('\nAre you satisfied with these kernel selections? (y/n): ', 's');
end

%% Define parameter ranges and sample parameter space
param_ranges.theta_cap = [1e-3, 1e-4, 1e-5];      % [hard, medium, easy]
param_ranges.kernel_size = [0.16, 0.04, 0.01];    % area ratio [hard, medium, easy]
param_ranges.SNR = [1.0, 3.16, 10.0];             % [hard, medium, easy]

[param_sets, descriptions] = sample_parameter_space(param_ranges);

%% Fixed parameters
fixed_params.num_kernels = num_kernels;
fixed_params.n = 1;  % number of energy layers per kernel
fixed_params.image_size = [500, 500];
fixed_params.relative_theta_cap = ones(1, num_kernels);  % [1, 1] for 2 kernels
fixed_params.relative_kernel_size = ones(1, num_kernels); % [1, 1] for 2 kernels

%% Generate datasets for each parameter set
num_sets = size(param_sets, 1);
datasets = struct('Y', {}, 'Y_clean', {}, 'A0', {}, 'A0_noiseless', {}, ...
                 'X0', {}, 'params', {}, 'b0', {});

%% Loop through parameter sets
for i = 1:num_sets
    confirmed = false;
    
    % Create figure once for this dataset
    fig = figure('Name', sprintf('Dataset %d/%d', i, num_sets), ...
                'Position', [100, 100, 1200, 800], ...
                'WindowState', 'normal');
    
    while ~confirmed
        fprintf('\nGenerating dataset %d/%d...\n', i, num_sets);
        
        % Extract parameters for this set
        theta_cap = param_sets(i,1);
        norm_kernel_size = param_sets(i,2);
        SNR = param_sets(i,3);
        
        % Calculate length ratio from area ratio
        length_ratio = sqrt(norm_kernel_size);
        
        % Calculate actual kernel size
        base_kernel_size = round(length_ratio * fixed_params.image_size(1));
        kernel_size = zeros(num_kernels, 2);
        for k = 1:num_kernels
            scaled_size = round(base_kernel_size * fixed_params.relative_kernel_size(k));
            kernel_size(k,:) = [scaled_size, scaled_size];
        end
        
        % Generate dataset
        [Y, Y_clean, A0, A0_noiseless, X0, b0] = generate_dataset_i(...
            theta_cap, norm_kernel_size, SNR, fixed_params, LDoS_sim, sliceidx);
        
        % Update display in existing figure
        clf(fig);  % Clear figure but keep it open
        imagesc(Y);
        title(sprintf('Dataset %d/%d\n𝜃_{cap}=%.0e, k_{size}=%.2f, SNR=%.1f', ...
              i, num_sets, theta_cap, norm_kernel_size, SNR));
        colorbar;
        axis square;
        
        % Center of screen, half size
        set(fig, 'Units', 'Normalized', ...
                'OuterPosition', [0.2, 0.2, 0.8, 0.8], ...
                'WindowState', 'normal');
        
        % Add buttons
        regenerate_btn = uicontrol('Style', 'pushbutton', ...
                                 'String', 'Regenerate', ...
                                 'Position', [20 20 100 40], ...
                                 'Callback', @(src,event) regenerate_dataset());
        
        confirm_btn = uicontrol('Style', 'pushbutton', ...
                               'String', 'Confirm', ...
                               'Position', [140 20 100 40], ...
                               'Callback', @(src,event) confirm_dataset());
        
        % Ensure figure is visible and properly sized
        drawnow;
        
        % Wait for button press
        uiwait(fig);
    end
    
    % Store results and close figure
    datasets(i) = store_dataset(Y, Y_clean, A0, A0_noiseless, X0, b0, theta_cap, kernel_size, SNR);
    fprintf('Dataset %d confirmed and saved.\n', i);
    close(fig);
end

% Helper function to store dataset
function dataset = store_dataset(Y, Y_clean, A0, A0_noiseless, X0, b0, theta_cap, kernel_size, SNR)
    dataset.Y = Y;
    dataset.Y_clean = Y_clean;
    dataset.A0 = A0;
    dataset.A0_noiseless = A0_noiseless;
    dataset.X0 = X0;
    dataset.params = struct('theta_cap', theta_cap, ...
                          'kernel_size', kernel_size, ...
                          'SNR', SNR);
    dataset.b0 = b0;
end

%% Save results
timestamp = datestr(now, 'yyyymmdd_HHMMSS');
save_filename = sprintf('synthetic_datasets_%s.mat', timestamp);
save(save_filename, 'datasets', 'param_sets', 'descriptions', ...
    'fixed_params', 'sliceidx', 'LDoS_sim');

fprintf('\nDatasets saved to: %s\n', save_filename); 

%% Display all generated datasets
fprintf('\nDisplaying all observations. Use scroll wheel to navigate.\n');
fprintf('Parameter settings for each slice:\n');

% Get number of datasets and create 3D array
num_sets = length(datasets);
[height, width] = size(datasets(1).Y);
all_observations = zeros(height, width, num_sets);

% Stack all observations and display parameters
for i = 1:num_sets
    all_observations(:,:,i) = datasets(i).Y;
    fprintf('Slice %d: 𝜃=%.0e, k=%.2f, SNR=%.1f\n', i, ...
        datasets(i).params.theta_cap, ...
        datasets(i).params.kernel_size(1), ...
        datasets(i).params.SNR);
end

% Display using d3gridDisplay
d3gridDisplay(all_observations, 'dynamic');
title('All Synthetic Observations');

%% Helper function to generate a single dataset
function [Y, Y_clean, A0, A0_noiseless, X0, b0] = generate_dataset_i(...
    theta_cap, norm_kernel_size, SNR, fixed_params, LDoS_sim, sliceidx)
    
    % Initialize kernels
    A0 = cell(1, fixed_params.num_kernels);
    A0_noiseless = cell(1, fixed_params.num_kernels);
    
    % Calculate length ratio from area ratio
    length_ratio = sqrt(norm_kernel_size);
    
    % Calculate average variance for universal noise
    avg_var = 0;
    for k = 1:fixed_params.num_kernels
        % Use length ratio for resizing
        kernel_length = round(length_ratio * fixed_params.image_size(1));
        A0_noiseless{k} = imresize(LDoS_sim(:,:,sliceidx(k)), [kernel_length, kernel_length]);
        A0_noiseless{k} = proj2oblique(A0_noiseless{k});
        avg_var = avg_var + var(A0_noiseless{k},0,"all");
    end
    avg_var = avg_var / fixed_params.num_kernels;
    
    % Apply universal noise
    eta_kernel = avg_var/SNR;
    for k = 1:fixed_params.num_kernels
        A0{k} = A0_noiseless{k} + sqrt(eta_kernel)*randn(size(A0_noiseless{k}));
        A0{k} = proj2oblique(A0{k});
    end
    
    % Generate activation maps
    theta = theta_cap/2 + theta_cap/2 * rand(1, fixed_params.num_kernels);
    X0 = zeros([fixed_params.image_size fixed_params.num_kernels]);
    
    for k = 1:fixed_params.num_kernels
        X0_good = false;
        while ~X0_good
            X0(:,:,k) = double(rand(fixed_params.image_size) <= theta(k));
            X0_good = sum(X0(:,:,k) ~= 0) > 0;
        end
    end
    
    % Generate observations
    b0 = randn;
    Y_clean = zeros(fixed_params.image_size);
    for k = 1:fixed_params.num_kernels
        Y_clean = Y_clean + convfft2(A0_noiseless{k}, X0(:,:,k));
    end
    
    % Add noise to match target SNR
    eta = var(Y_clean, 0, "all") / SNR;
    Y = Y_clean + sqrt(eta) * randn(size(Y_clean));
end

% Callback function for regenerate button
function regenerate_dataset()
    uiresume;
end

% Callback function for confirm button
function confirm_dataset()
    assignin('base', 'confirmed', true);
    uiresume;
end

