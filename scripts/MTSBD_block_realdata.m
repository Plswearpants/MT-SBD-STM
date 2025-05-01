%% Block 1: Load the .3ds data

% INPUTS
% 1: Data file to load, including file type ('QPI.3ds' for example)
% 2: Smoothing sigma for current data

% OUTPUTS
% header: Variable containing all experimental parameters
% I: Current data, smoothed by sigma
% dIdV: Numerically differentiated current data
% voltage: Vector of voltages for current
% midV: Vector on voltages for dIdV/QPI (midpoint of voltage vector)
% QPI: Fourier transformed dIdV data

% Modified function load3dsall from supplied matlab code from Nanonis
[header, par, I, dIdV, LockindIdV, bias, midV, QPI, LockinQPI] = load3dsall('Grid Spectroscopy_QPI_24nm_3days002.3ds', 5);
xsize = header.grid_dim(1);
ysize = header.grid_dim(2);
elayer = header.points;
estart = par(1);
eend = par(2);
energy_range = linspace(estart, eend, elayer);

Y = dIdV;
num_slices = size(Y,3);
spatial = size(Y,1);
%% Block 2: Data preprocessing

%% 2.1: Remove bragg peaks
[Y]=removeBragg(Y);

%% 2.2: crop dataset
mask= maskSquare(Y,0,301);
Y= gridCropMask(Y, mask);

%% 2.3 defect masking

f1=figure;
d3gridDisplay(Y,'dynamic');
index = input('Enter defect slice number: ');
close(f1);

% methods: 
% 1. Gaussian window "gw"
% 2. truncated gaussian gaussian smoothing "tg"
% 3. thresholding and remove defect features "threshold"

method = 'tg';

switch method
    case 'gw'
        % Apply Gaussian window masking
        [Y_masked, ~] = defect_masking(Y, index);
    case 'tg'
        % Apply flat disk mask with Gaussian smoothing
        [Y_masked, defect_mask] = gaussianMaskDefects(Y,index);
    case 'threshold'
        % Apply threshold-based defect masking
        [Y_masked, defect_mask] = thresholdDefects(Y, index);
    otherwise
        error('Unknown defect masking method. Choose "gw", "disk", or "threshold".');
end
%% Block 3 data selection 
rangetype='dynamic';
f1=figure;
d3gridDisplay(Y,rangetype);
params.start_slice = input('Enter QPI starting slice number: ');
params.end_slice = input('Enter QPI ending slice number: ');
Y=Y(:,:,params.start_slice:params.end_slice);
close(f1);
% update energy range from start and end slice
energy_range = energy_range(params.start_slice:params.end_slice);


%%
f2=figure;
d3gridDisplay(Y,rangetype);
params.ref_slice = input('Enter reference slice number: ');
num_kernels = input('enter the number of kernels you wish to apply: ');
close(f2);


% Validate the input
if isempty(params.ref_slice) || ~isnumeric(params.ref_slice) || params.ref_slice < 1 || params.ref_slice > size(Y, 3)
    error('Invalid reference slice. Please enter a number between 1 and %d.', size(Y, 3));
end
fprintf('Using slice %d as reference\n', params.ref_slice);

Y = normalizeBackgroundToZeroMean3D(Y, 'dynamic', params.ref_slice);  % normalize Y
Y = proj2oblique(Y);

% Extract reference slice and initialize kernel
Y_ref = Y(:,:,params.ref_slice);
Y_ref = normalizeBackgroundToZeroMean3D(Y_ref, 'dynamic');  % normalize Y
Y_ref = proj2oblique(Y_ref);

% Display the selected reference slice
figure;
imagesc(Y_ref);
colorbar;
title(sprintf('Reference Slice %d', params.ref_slice));
axis square;

% Initialize kernels
% draw square on the data to include as many visible ripples of the scattering as possible 
same_size = 1;
kerneltype = 'selected';
window_type = {'gaussian', 1};
%window_type = '';

if same_size
    [square_size] = squareDrawSize(Y_ref);
    kernel_sizes = repmat(square_size,[num_kernels,1]);
    A1 = initialize_kernels(Y_ref, num_kernels, kernel_sizes, kerneltype, window_type);
else
    A1 = cell(1, num_kernels);
    kernel_sizes = zeros(num_kernels, 2); % Store sizes of each kernel [height, width]
    for k = 1:num_kernels
        fprintf('Select region for kernel %d/%d\n', k, num_kernels);
        [square_size,position, mask] = squareDrawSize(target_data);           	% determine kernel size
        [A1{k}, ~] = gridCropMask(target_data, mask);   % the cropped real data as kernel
        % Need to put each slice back onto the sphere
        A1{k} = proj2oblique(A1{k});
        % Store the kernel size
        kernel_sizes(k,:) = size(A1{k});
    end
end

%% (ESS)noise level determination 
eta_data = estimate_noise(Y_ref, 'std');  

%% (Opt) determine SNR
SNR_data= var(A1{1,1}(:))/eta_data;
fprintf('SNR_data = %d', SNR_data);

%% Block 2: Find Optimal Activation for Reference Slice
% Set up display functions
figure;
dispfun = cell(1, num_kernels);
for n = 1:num_kernels
    dispfun{n} = @(Y, A, X, kernel_sizes, kplus) showims(Y_ref, A1{n}, X, A, X, kernel_sizes, kplus, 1);
end

% SBD settings.
miniloop_iteration = 2;
outerloop_maxIT= 20;

params.lambda1 = [1e-1, 1e-1,0.1,0.1];  % regularization parameter for Phase I
%params.lambda1 = [0.15, 0.15, 0.15, 0.15, 0.15];  % regularization parameter for Phase I
params.phase2 = false;
params.kplus = ceil(0.5 * kernel_sizes);
params.lambda2 = [2e-2, 2e-2];  % FINAL reg. param. value for Phase II
params.nrefine = 3;
params.signflip = 0.2;
params.xpos = true;
params.getbias = true;
params.Xsolve = 'FISTA';

% noise variance for computeResidualQuality.m
params.noise_var = eta_data;

%% Run and save 
% 2. The fun part
[Aout, Xout, bout, extras] = MT_SBD(Y_ref, kernel_sizes, params, dispfun, A1, miniloop_iteration, outerloop_maxIT);
%% 
% Choose which demixing method to use
params.use_Xregulated = false;  % Set to true to use MTSBD_Xregulated, false for SBD_test_multi_demixing

% Run SBD on reference slice
fprintf('Finding optimal activation for reference slice %d...\n', params.ref_slice);
if params.use_Xregulated
    [REG_A_ref, REG_X_ref, REG_bout, REG_extras] = MTSBD_Xregulated(...
        Y_ref, kernel_sizes, params, dispfun, A1, initial_iteration, maxIT);
    % Store reference results
    REG_ref_results = struct();
    REG_ref_results.A = REG_A_ref;
    REG_ref_results.X = REG_X_ref;
    REG_ref_results.extras = REG_extras;
else
    [A_ref, X_ref, bout, extras] = SBD_test_multi_demixing(...
        Y_ref, kernel_sizes, params, dispfun, A1, initial_iteration, maxIT);
    % Store reference results
    ref_results = struct();
    ref_results.A = A_ref;
    ref_results.X = X_ref;
    ref_results.extras = extras;
end



%% Block 3: Find Most Isolated Points and Initialize Kernels
fprintf('Calculating isolation scores and finding most isolated points...\n');

% Choose target kernel sizes first
type = 'kernel_sizes_cap';
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

params.use_Xregulated = true;


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
    [REG_Aout_ALL, REG_Xout_ALL, REG_bout_ALL, REG_extras_ALL] = MTSBD_Xregulated_all_slices(...
        Y_used, kernel_sizes_used, params, dispfun, A1_used, initial_iteration, maxIT);
else
    [Aout_slice, Xout_slice, bout_slice, slice_extras] = MTSBD_demixing_all_slice(...
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
