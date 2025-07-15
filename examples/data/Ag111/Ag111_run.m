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
[header, par, I, dIdV, LockindIdV, bias, midV, QPI, LockinQPI] = load3dsall('Grid Spectroscopy002.3ds', 5);
xsize = header.grid_dim(1);
ysize = header.grid_dim(2);
elayer = header.points;
estart = par(1);
eend = par(2);
energy_range = linspace(estart, eend, elayer);
data_original = dIdV;
num_slices = size(data_original,3);
spatial = size(data_original,1);

%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Block 2: Data preprocessing
data_carried = data_original;
rangetype='dynamic';
figure;
d3gridDisplay(data_carried,rangetype);
slice_normalize = input('slice to normalize: ');
%% 2.2b defect masking
% Normalize background 
[data_carried] = normalizeBackgroundToZeroMean3D(data_carried, rangetype, slice_normalize);

f1=figure;
d3gridDisplay(data_carried,'dynamic');
index = input('Enter defect slice number: ');
num_defect_type = input('enter how many types of defects to mask: ');
close(f1);
% methods: 
% 1. Gaussian window "gw"
% 2. truncated gaussian gaussian smoothing "tg"
% 3. thresholding and remove defect features "threshold"
method = 'tg';

switch method
    case 'gw'
        % Apply Gaussian window masking
        [data_masked, ~] = defect_masking(data_carried, index);
    case 'tg'
        % Apply flat disk mask with Gaussian smoothing
        [data_masked, defect_loc] = gaussianMaskDefects(Y,index, num_defect_type);
    case 'threshold'
        % Apply threshold-based defect masking
        [data_masked, defect_mask] = thresholdDefects(data_carried, index);
    otherwise
        error('Unknown defect masking method. Choose "gw", "disk", or "threshold".');
end
data_carried = data_masked;

%% 2.3a: Correct streak
% Normalize background 
[data_carried] = normalizeBackgroundToZeroMean3D(data_carried, rangetype, slice_normalize);

[data_streakremoved, QPI_nostreaks] = RemoveStreaks(data_carried, 'Direction', 'vertical');
data_carried = data_streakremoved;

%% 2.4: crop dataset
mask= maskSquare(data_carried,0,slice_normalize);
data_cropped= gridCropMask(data_carried, mask);
data_carried = data_cropped;

%% 2.end: Normalize background 
[Y] = normalizeBackgroundToZeroMean3D(data_carried, rangetype, slice_normalize);

%% 3pre: Save the preprocessed data~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
save('Ag002_preprocessed_selected.mat', 'Y', 'data_masked', 'data_cropped', 'data_streakremoved', 'data_original', "energy_range",'energy_selected','defect_loc')


%% Block 3 data selection 
rangetype='dynamic';
f1=figure;
d3gridDisplay(Y,rangetype);
params.slices = input('input a list of slices: ');
Y=data_carried(:,:,params.slices);
energy_selected = energy_range(params.slices);
close(f1);
num_slices = size(Y,3);

%% 3.end: Normalize background 
[Y] = normalizeBackgroundToZeroMean3D(Y, rangetype, 1);

%% 4_pre Pick reference slice and Initialize reference kernels
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

%% Initialize reference kernels
% draw square on the data to include as many visible ripples of the scattering as possible 
same_size = 1;
kerneltype = 'selected';
window_type = {'gaussian', 2.5};
%window_type = '';


if same_size
    [square_size] = squareDrawSize(Y_ref);
    %square_size=[80,80];
    kernel_sizes = repmat(square_size,[num_kernels,1]);
    A1_ref = initialize_kernels(Y_ref, num_kernels, kernel_sizes, kerneltype, window_type);
else
    A1_ref = cell(1, num_kernels);
    kernel_sizes = zeros(num_kernels, 2); % Store sizes of each kernel [height, width]
    for k = 1:num_kernels
        fprintf('Select region for kernel %d/%d\n', k, num_kernels);
        [square_size,position, mask] = squareDrawSize(Y_ref);           	% determine kernel size
        [A1_ref{k}, ~] = gridCropMask(Y_ref, mask);   % the cropped real data as kernel
        % Need to put each slice back onto the sphere
        A1_ref{k} = proj2oblique(A1_ref{k});
        % Store the kernel size
        kernel_sizes(k,:) = size(A1_ref{k});
    end
end

%% (ESS)noise level determination 
eta_data = estimate_noise(Y_ref, 'std');  

%% Block 4: Find Optimal Activation for Reference Slice
% Set up display functions
figure;
dispfun = cell(1,  num_kernels);
for n = 1:num_kernels
    dispfun{n} = @(Y, A, X, kernel_sizes, kplus) showims(Y_ref, A1_ref{n}, X, A, X, kernel_sizes, kplus, 1);
end

% SBD settings.
miniloop_iteration = 3;
outerloop_maxIT= 15;
%params_ref.energy = energy_selected(params.ref_slice);
params_ref.lambda1 = [0.05, 0.03,0.02, 0.05];  % regularization parameter for Phase I
%params_ref.lambda1 = [0.15, 0.15, 0.15, 0.15, 0.15];  % regularization parameter for Phase I
params_ref.phase2 = true;
params_ref.kplus = ceil(0.2 * kernel_sizes);
params_ref.lambda2 = [0.05, 0.05, 0.05, 0.15];  % FINAL reg. param. value for Phase II
params_ref.nrefine = 4;
params_ref.signflip = 0.2;
params_ref.xpos = true;
params_ref.getbias = true;
params_ref.Xsolve = 'FISTA';

% noise variance for computeResidualQuality.m
params_ref.noise_var = eta_data;

% Run and save 
% 2. The fun part
[A_ref, X_ref, b_ref, extras_ref] = MT_SBD(Y_ref, kernel_sizes, params_ref, dispfun, A1_ref, miniloop_iteration, outerloop_maxIT);

%% Visualize Reference result 
visualizeRealResult(Y_ref,A_ref, X_ref, b_ref, extras_ref);

%% Pad the A_ref to be size defined by user, normalize and use them as the A1
target_size = [100, 100];
kernel_sizes_pad = repmat(target_size,[num_kernels,1]);
%kernel_sizes_pad = [[120,120];[120,120];[65,65]];
A_pre_pad = A_ref;
A1_ref = cell(1, num_kernels);
for k = 1:num_kernels
    sz = size(A_pre_pad{k});
    pad_h = kernel_sizes_pad(k,1) - sz(1);
    pad_w = kernel_sizes_pad(k,2) - sz(2);

    % Calculate pre- and post-padding for centering
    pre_h = floor(pad_h / 2);
    post_h = ceil(pad_h / 2);
    pre_w = floor(pad_w / 2);
    post_w = ceil(pad_w / 2);

    % Pad so that the kernel is centered
    A1_ref{k} = padarray(A_pre_pad{k}, [pre_h, pre_w], 'pre');
    A1_ref{k} = padarray(A1_ref{k}, [post_h, post_w], 'post');
    
    A1_ref{k} = proj2oblique(A1_ref{k});
end

% visualize the padded kernels
figure;
for k = 1:num_kernels
    subplot(1,num_kernels,k);
    imagesc(A1_ref{k}); axis square;
    colorbar;
end

% define the intial activation map using X_ref
for k = 1:num_kernels
    params_ref.xinit{k}.X = X_ref(:,:,k);
    params_ref.xinit{k}.b = extras_ref.phase1.biter(k); 
end
%% Set ups before padded run 
% Set up display functions
figure;
dispfun = cell(1, num_kernels);
for n = 1:num_kernels
    dispfun{n} = @(Y, A, X, kernel_sizes, kplus) showims(Y_ref, A1_ref{n}, X, A, X, kernel_sizes, kplus, 1);
end

% SBD settings.
miniloop_iteration = 1;
outerloop_maxIT= 5;

params_ref.lambda1 = [0.03, 0.03, 0.03,0.03];  % regularization parameter for Phase I
%params_ref.lambda1 = [0.15, 0.15, 0.15, 0.15, 0.15];  % regularization parameter for Phase I
params_ref.phase2 = false;
params_ref.kplus = ceil(0.5 * kernel_sizes);
params_ref.lambda2 = [2e-2, 2e-2];  % FINAL reg. param. value for Phase II
params_ref.nrefine = 3;
params_ref.signflip = 0.2;
params_ref.xpos = true;
params_ref.getbias = true;
params_ref.Xsolve = 'FISTA';

% noise variance for computeResidualQuality.m
params_ref.noise_var = eta_data;
%% Run the padded initialization 
% 2. The fun part
[A_ref_pad, X_ref_pad, b_ref_pad, extras_ref_pad] = MT_SBD(Y_ref, kernel_sizes_pad, params_ref, dispfun, A1_ref, miniloop_iteration, outerloop_maxIT);

%% Visualize Padded result 
visualizeRealResult(Y_ref,A_ref_pad, X_ref_pad, b_ref_pad, extras_ref_pad);

%% Block 3: Find Most Isolated Points and Initialize Kernels
num_slices = size(Y,3);

% Choose method for kernel center selection
fprintf('Choose method for kernel center selection:\n');
fprintf('1. Find most isolated points automatically\n');
fprintf('2. Manually select 3 kernel centers\n');
choice = input('Enter choice (1 or 2): ');

if choice == 1
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
        %threshold = max(X_ref(:,:,k),[],'all')/10;
        threshold = 0;
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

elseif choice == 2
    fprintf('Manual kernel center selection mode...\n');
    
    % Display the reference image for manual selection
    figure('Name', 'Manual Kernel Center Selection');
    imagesc(Y_ref);
    title('Click on 3 points to select kernel centers. Press Enter when done.');
    colormap(gray);
    colorbar;
    
    % Get user clicks for kernel centers
    kernel_centers = zeros(num_kernels, 2);
    for k = 1:num_kernels
        fprintf('Click on center for kernel %d/%d\n', k, num_kernels);
        [x, y] = ginput(1);
        kernel_centers(k,:) = [round(y), round(x)]; % Convert to [row, col] format
        
        % Mark the selected point
        hold on;
        scatter(x, y, 100, 'r', '*');
        text(x+5, y+5, sprintf('K%d', k), 'Color', 'red', 'FontSize', 12, 'FontWeight', 'bold');
        hold off;
    end
    
    fprintf('Kernel centers selected:\n');
    for k = 1:num_kernels
        fprintf('Kernel %d: (%d, %d)\n', k, kernel_centers(k,1), kernel_centers(k,2));
    end
    
    % Use target kernel sizes from the reference slice
    target_kernel_sizes = kernel_sizes;
    
else
    error('Invalid choice. Please enter 1 or 2.');
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
    imagesc(A1_all{params.ref_slice,k});
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

%% Convert A1_all to matrix form and prepare noise
A1_all_matrix = cell(num_kernels,1);
for k = 1:num_kernels
    A1_all_matrix{k} = zeros(size(A1_all{1,k},1),size(A1_all{1,k},2),num_slices);
    for s = 1:num_slices
        A1_all_matrix{k}(:,:,s) = A1_all{s,k};
    end
end

eta_data3d = estimate_noise3D(Y, 'std');  

%% initialize xinit for all slices with reference slice
for k = 1:num_kernels
    params.xinit{k}.X = X_ref(:,:,k);
    b_temp = extras_ref.phase1.biter(k); 
    params.xinit{k}.b = repmat(b_temp,[num_slices,1]);
end


%% Run for all_slice
% SBD settings.
miniloop_iteration = 1;
outerloop_maxIT= 8;

params.lambda1 = [0.05, 0.03, 0.03];  % regularization parameter for Phase I
%params.lambda1 = [0.15, 0.15, 0.15, 0.15, 0.15];  % regularization parameter for Phase I
params.phase2 = false;
params.kplus = ceil(0.5 * kernel_sizes);
params.lambda2 = [2e-2, 2e-2];  % FINAL reg. param. value for Phase II
params.nrefine = 3;
params.signflip = 0.2;
params.xpos = true;
params.getbias = true;
params.Xsolve = 'FISTA';
params.use_Xregulated = false;
params.noise_var = eta_data3d;

kernel_sizes_single = squeeze(max(kernel_sizes,[],1));
Y_used = Y;
A1_used = A1_all_matrix;

% Set up display functions
figure;
dispfun = cell(1, num_kernels);
for n = 1:num_kernels
    dispfun{n} = @(Y, A, X, kernel_sizes_sing, kplus) showims(Y_used, A1_used{n}, X_ref(:,:,n), A, X, kernel_sizes_single, kplus, 1);
end

% Update params for MTSBD
for k = 1:num_kernels
    params.xinit{k}.X = X_ref(:,:,k);
    b_temp = extras_ref.phase1.biter(k); 
    params.xinit{k}.b = repmat(b_temp,[num_slices,1]);
end

if params.use_Xregulated
    [REG_Aout_ALL, REG_Xout_ALL, REG_bout_ALL, REG_extras_ALL] = MTSBD_Xregulated_all_slices(...
        Y_used, kernel_sizes, params, dispfun, A1_used, miniloop_iteration, outerloop_maxIT);
else
    [Aout_ALL, Xout_ALL, bout_ALL, ALL_extras] = MTSBD_all_slice(...
        Y_used, kernel_sizes, params, dispfun, A1_used, miniloop_iteration, outerloop_maxIT);
end
%%
save('Ag111_slices.mat', 'Y', 'A1_all_matrix','Aout_ALL', 'Xout_ALL', 'bout_ALL', 'ALL_extras');