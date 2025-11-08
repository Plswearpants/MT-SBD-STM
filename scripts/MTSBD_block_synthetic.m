%% Block 1: Generate Test Data and Initialize Kernels
% Generate test set
SNR = 2;              % Signal-to-noise ratio
N_obs = 100;           % Observation lattice size 
observation_resolution = 3;  % Resolution: 3 pixels per lattice site
defect_density = 0.005;    
num_slices = 2;

% Generate data with custom parameters
[Y, A0_noiseless, X0, params] = properGen_full(SNR, N_obs, observation_resolution, defect_density, ...
    'LDoS_path', 'example_data/LDoS_single_defects_self=0.6_save.mat', ...
    'num_slices', num_slices, ...
    'vis', 0);

% Extract parameters from generated data
A0 = params.A0;
num_kernels = params.num_kernels;
num_slices = params.num_slices; 

% Standardize kernel sizes structure
kernel_sizes = zeros(num_slices, num_kernels, 2);
for s = 1:num_slices
    for k = 1:num_kernels
        kernel_sizes(s,k,:) = size(A0{s,k});
    end
end

% Normalize and visualize
Y = normalizeBackgroundToZeroMean3D(Y, 'dynamic'); 
Y = proj2oblique(Y);

% Choose reference slice selection method

% Display the 3D data and let the user input the reference slice
figure; 
d3gridDisplay(Y,'dynamic');
title('Select a reference slice');
params.ref_slice = input('Enter reference slice number: ');

% Validate the input
if isempty(params.ref_slice) || ~isnumeric(params.ref_slice) || params.ref_slice < 1 || params.ref_slice > size(Y, 3)
    error('Invalid reference slice. Please enter a number between 1 and %d.', size(Y, 3));
end
fprintf('Using slice %d as reference\n', params.ref_slice);

% Extract reference slice data
Y_ref = Y(:,:,params.ref_slice);
X0_ref = X0;
A0_ref = cell(1, num_kernels);
for i = 1:num_kernels
    A0_ref{i} = A0{params.ref_slice,i};
end

% Display the selected reference slice
figure;
imagesc(Y_ref);
colorbar;
title(sprintf('Reference Slice %d', params.ref_slice));
axis square;

%% Initialize kernels
window_type = {'gaussian', 2.5};  % Example: gaussian window with alpha
kernel_sizes_ref = reshape(kernel_sizes(params.ref_slice,:,:),[num_kernels,2]);
figure;
for k = 1: num_kernels 
    subplot(1,num_kernels,k); imagesc(A0{params.ref_slice,k})
    axis square
end
%
A1 = initialize_kernels(Y_ref, num_kernels, kernel_sizes_ref, 'selected', window_type);

% Display initialized kernels
figure;
for n = 1:num_kernels
    subplot(1, num_kernels, n);
    imagesc(A1{n});
    title(sprintf('Initial Kernel %d', n));
    colorbar;
    axis square;
end

%% Block 2: Find Optimal Activation for Reference Slice
% Set up display functions
figure;
dispfun = cell(1, num_kernels);
for n = 1:num_kernels
    dispfun{n} = @(Y, A, X, kernel_sizes, kplus) showims(Y_ref, A1{n}, X0_ref(:,:,n), A, X, kernel_sizes, kplus, 1);
end

% Choose which demixing method to use
params.use_Xregulated = false;  % Set to true to use MTSBD_Xregulated, false for SBD_test_multi_demixing

% SBD settings
initial_iteration = 1;
maxIT = 15;
params.phase2 = false;  % whether Enable Phase 2
params.lambda1 = [3e-2, 3e-2, 3e-2, 3e-2];  % regularization parameter for Phase I
params.kplus = ceil(0.5 * kernel_sizes_ref);
params.lambda2 = [1e-2, 1e-2, 1e-2, 1e-2];  % FINAL reg. param. value for Phase II
params.nrefine = 5;  % Number of refinements for Phase II
params.signflip = 0.2;
params.xpos = true;
params.getbias = true;
params.Xsolve = 'FISTA';
params.X0 = X0_ref;
params.A0 = A0_ref;
params.xinit = [];
% Add cross-correlation regularization parameter
params.gamma = 5e-2;  % Cross-correlation regularization parameter

% Run SBD on reference slice
fprintf('Finding optimal activation for reference slice %d...\n', params.ref_slice);
if params.use_Xregulated
    [REG_A_ref, REG_X_ref, REG_bout, REG_extras] = MTSBD_synthetic_Xregulated(...
        Y_ref, kernel_sizes_ref, params, dispfun, A1, initial_iteration, maxIT);
    % Store reference results
    REG_ref_results = struct();
    REG_ref_results.A = REG_A_ref;
    REG_ref_results.X = REG_X_ref;
    REG_ref_results.extras = REG_extras;
else
    [A_ref, X_ref, bout, extras] = MTSBD_synthetic(...
        Y_ref, kernel_sizes_ref, params, dispfun, A1, initial_iteration, maxIT);
    % Store reference results
    ref_results = struct();
    ref_results.A = A_ref;
    ref_results.X = X_ref;
    ref_results.extras = extras;
end

%% Block 3: Find Most Isolated Points and Initialize Kernels
fprintf('Calculating isolation scores and finding most isolated points...\n');

% Choose target kernel sizes first
type = 'kernel_sizes_all';
if strcmp(type, 'ref_kernel_sizes')
    target_kernel_sizes = squeeze(kernel_sizes(params.ref_slice,:,:));
elseif strcmp(type, 'kernel_sizes_cap')
    target_kernel_sizes = squeeze(max(kernel_sizes,[],1));
elseif strcmp(type, 'kernel_sizes_all')
    target_kernel_sizes = kernel_sizes;
end

% Get defect positions from X_ref
most_isolated_points = cell(1, num_kernels);
isolation_scores = cell(1, num_kernels);
defect_positions = cell(1, num_kernels);
num_defects = zeros(1, num_kernels);

% Create figure for activation value distributions
figure('Name', 'Activation Value Distributions');
for k = 1:num_kernels
    % Plot histogram of activation values
    subplot(2, num_kernels, k);
    activation_values = X_ref(:,:,k);
    h = histogram(activation_values(activation_values > 0), 50);
    set(gca, 'YScale', 'log');
    title(sprintf('Kernel %d Activation Distribution', k));
    xlabel('Activation Value');
    ylabel('Frequency (log scale)');
    
    % Add vertical line for threshold
    threshold = max(X_ref(:,:,k),[],'all')/10;
    hold on;
    xline(threshold, 'r--', 'Threshold');
    hold off;
    
    % Plot cumulative distribution
    subplot(2, num_kernels, k + num_kernels);
    [counts, edges] = histcounts(activation_values(activation_values > 0), 50, 'Normalization', 'cdf');
    stairs(edges(1:end-1), counts);
    title(sprintf('Kernel %d Cumulative Distribution', k));
    xlabel('Activation Value');
    ylabel('Cumulative Frequency');
    hold on;
    xline(threshold, 'r--', 'Threshold');
    hold off;
    
    % Get positions of defects above threshold
    [rows, cols] = find(X_ref(:,:,k) > threshold);
    defect_positions{k} = [rows, cols];
    num_defects(k) = size(defect_positions{k}, 1);
    fprintf('Kernel %d has %d defects above threshold %.4f\n', k, num_defects(k), threshold);
end

% Calculate isolation scores for each kernel
for k = 1:num_kernels
    if num_defects(k) == 0
        warning('No defects found for kernel %d', k);
        continue;
    end
    
    % Create summed activation map of all other kernels
    X_others = zeros(size(X_ref(:,:,1)));
    for l = 1:num_kernels
        if l ~= k
            X_others = X_others + X_ref(:,:,l);
        end
    end
    
    % Get positions of defects in other kernels
    [other_rows, other_cols] = find(X_others > max(X_others,[],'all')/10);
    other_positions = [other_rows, other_cols];
    
    % Use target kernel size for boundary check
    half_kernel_size = floor(target_kernel_sizes(k,:)/2);
    
    % First filter out points too close to boundaries
    valid_points = true(num_defects(k), 1);
    for i = 1:num_defects(k)
        y = defect_positions{k}(i,1);
        x = defect_positions{k}(i,2);
        
        if y <= half_kernel_size(1) || y >= size(X_ref,1) - half_kernel_size(1) || ...
           x <= half_kernel_size(2) || x >= size(X_ref,2) - half_kernel_size(2)
            valid_points(i) = false;
        end
    end
    
    % Only proceed with valid points
    valid_defects = defect_positions{k}(valid_points,:);
    if isempty(valid_defects)
        error('No valid isolated points found for kernel %d - all points are too close to image boundaries', k);
    end
    
    % Calculate isolation scores only for valid points
    S_k = zeros(size(valid_defects, 1), 1);
    for i = 1:size(valid_defects, 1)
        diffs = other_positions - valid_defects(i,:);
        distances = sum(diffs.^2, 2);
        S_k(i) = min(distances);
    end
    
    % Find the most isolated point among valid points
    [max_score, max_idx] = max(S_k);
    valid_indices = find(valid_points);
    max_idx = valid_indices(max_idx);
    
    most_isolated_points{k} = defect_positions{k}(max_idx,:);
    isolation_scores{k} = S_k;
    
    fprintf('Kernel %d: Most isolated point at (%d,%d) with score %.2f\n', ...
        k, most_isolated_points{k}(1), most_isolated_points{k}(2), max_score);
end

% Visualize the results
figure('Name', 'Isolation Analysis');
for k = 1:num_kernels
    subplot(2, num_kernels, k);
    imagesc(X_ref(:,:,k));
    hold on;
    scatter(defect_positions{k}(:,2), defect_positions{k}(:,1), 50, 'w', 'o');
    if ~isempty(most_isolated_points{k})
        scatter(most_isolated_points{k}(2), most_isolated_points{k}(1), 100, 'r', '*');
    end
    title(sprintf('Kernel %d', k));
    colorbar;
    axis square;
    hold off;
    
    subplot(2, num_kernels, k + num_kernels);
    X_others = zeros(size(X_ref(:,:,1)));
    for l = 1:num_kernels
        if l ~= k
            X_others = X_others + X_ref(:,:,l);
        end
    end
    imagesc(X_others);
    hold on;
    if ~isempty(most_isolated_points{k})
        scatter(most_isolated_points{k}(2), most_isolated_points{k}(1), 100, 'r', '*');
    end
    title(sprintf('Other Kernels (not %d)', k));
    colorbar;
    axis square;
    hold off;
end

% Store isolation analysis results
isolation_analysis = struct();
isolation_analysis.defect_positions = defect_positions;
isolation_analysis.most_isolated_points = most_isolated_points;
isolation_analysis.isolation_scores = isolation_scores;
isolation_analysis.num_defects = num_defects;

% Prepare ground truth kernels
if strcmp(type, 'kernel_sizes_all')
    A0_used = A0;
else
    A0_used = padKernels(A0_noiseless, SNR, target_kernel_sizes);
end

% Align most isolated points with ground truth activation
[~, offset, ~] = alignActivationMaps(X0_ref, X_ref, kernel_sizes_ref);
used_most_isolated_points = cell(1, num_kernels);
for k = 1:num_kernels
    used_most_isolated_points{k} = most_isolated_points{k} + offset(k,:);
end

% Display most isolated points
figure;
imagesc(Y_ref);
hold on;
for k = 1:num_kernels
    scatter(most_isolated_points{k}(2), most_isolated_points{k}(1), 100, 'r', '*');
end

% Convert most_isolated_points to matrix format
kernel_centers = zeros(num_kernels, 2);
for k = 1:num_kernels
    if ~isempty(most_isolated_points{k})
        kernel_centers(k,:) = most_isolated_points{k};
    else
        error('No isolated point found for kernel %d', k);
    end
end

% Initialize kernels for all slices
A1_all = cell(num_slices, num_kernels);
matrix = true;  % Use matrix form for A1
change_size = false;

for s = 1:num_slices
    fprintf('Initializing kernels for slice %d/%d...\n', s, num_slices);
    if matrix
        A1_all(s,:) = initialize_kernels_proliferation(Y(:,:,s), num_kernels, kernel_centers, window_type, target_kernel_sizes, 'interactive', change_size);
    else 
        A1_all(s,:) = initialize_kernels_proliferation(Y(:,:,s), num_kernels, kernel_centers, window_type, squeeze(kernel_sizes(s,:,:)), 'interactive', change_size);
    end
end

% Visualize initialized kernels for reference slice
figure('Name', 'Initialized Kernels (Reference Slice)');
for k = 1:num_kernels
    subplot(1,num_kernels,k);
    imagesc(A1_all{params.ref_slice+1,k});
    colormap(gray);
    colorbar;
    title(sprintf('Initialized Kernel %d', k));
    axis square;
end

% Store initialization results
init_results = struct();
init_results.A = A1_all;
init_results.kernel_centers = kernel_centers;
init_results.kernel_sizes = target_kernel_sizes;

fprintf('Kernel initialization complete for all slices.\n');

%% Convert A1_all to matrix form
A1_all_matrix = cell(num_kernels,1);
for k = 1:num_kernels
    A1_all_matrix{k} = zeros(size(A1_all{1,k},1),size(A1_all{1,k},2),num_slices);
    for s = 1:num_slices
        A1_all_matrix{k}(:,:,s) = A1_all{s,k};
    end
end

%% Run for all_slice

params.use_Xregulated = false;


kernel_sizes_used = squeeze(max(kernel_sizes,[],1));
Y_used = Y;
A1_used = A1_all_matrix;

% Set up display functions
figure;
dispfun = cell(1, num_kernels);
for n = 1:num_kernels
    dispfun{n} = @(Y, A, X, kernel_sizes, kplus) showims(Y_used, A1_used{n}, X0_ref(:,:,n), A, X, kernel_sizes, kplus, 1);
end

% Update params for MTSBD
params.X0 = X0;
params.A0 = A0_used;
for k = 1:num_kernels
    params.xinit{k}.X = X_ref(:,:,k);
    b_temp = extras.phase1.biter(k); 
    params.xinit{k}.b = repmat(b_temp,[num_slices,1]);
end

if params.use_Xregulated
    [REG_Aout_ALL, REG_Xout_ALL, REG_bout_ALL, REG_extras_ALL] = MTSBD_synthetic_Xregulated_all_slices(...
        Y_used, kernel_sizes_used, params, dispfun, A1_used, initial_iteration, maxIT);
else
    [Aout_slice, Xout_slice, bout_slice, slice_extras] = MTSBD_synthetic_all_slice(...
        Y_used, kernel_sizes_used, params, dispfun, A1_used, initial_iteration, maxIT);
end

%% Visualization
% Convert Aout_slice to cell format for visualization
Aout_cell = cell(num_slices, num_kernels);
for s = 1:num_slices
    for k = 1:num_kernels
        Aout_cell{s,k} = Aout_slice{k}(:,:,s);
    end
end

% Visualize results for each slice
for s = 1:num_slices
    % Prepare slice-specific data
    Y_slice = Y(:,:,s);
    A0_slice = cell(1, num_kernels);
    for k = 1:num_kernels
        A0_slice{k} = A0_used{k}(:,:,s);
    end
    Aout_slice_cell = cell(1, num_kernels);
    for k = 1:num_kernels
        Aout_slice_cell{k} = Aout_cell{s,k};
    end
    
    % Call visualizeResults for this slice
    visualizeResults(Y_slice, A0_slice, Aout_slice_cell, X0, Xout_slice, bout_slice(s,:), slice_extras, [s, 1]);
end
