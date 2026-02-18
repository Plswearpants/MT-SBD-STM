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
[header, par, I, dIdV, LockindIdV, bias, midV, QPI, LockinQPI] = load3dsall('QPImap012.3ds', 10);
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
% initialize the preprocessing parameters
preprocessing_params = struct();
data_carried = data_original;
rangetype='dynamic';
figure;
d3gridDisplay(data_carried,rangetype);
preprocessing_params.slice_normalize = input('slice to normalize: ');

%% 2.1: Remove bragg peaks
% Normalize background 
[data_carried] = normalizeBackgroundToZeroMean3D(data_carried, rangetype, preprocessing_params.slice_normalize);
% Bragg remove
[data_braggremoved]=removeBragg(data_carried);
data_carried = data_braggremoved;

%% 2.2: crop dataset
mask= maskSquare(data_carried,0,40,'square');
data_cropped= gridCropMask(data_carried, mask);
data_carried = data_cropped;

%% 2.3 data selection 
rangetype='dynamic';
figure;
d3gridDisplay(data_carried,rangetype);
preprocessing_params.slices = input('input a list of slices: ');
data_selected=data_carried(:,:,preprocessing_params.slices);
energy_selected = energy_range(preprocessing_params.slices);
close;
data_carried = data_selected; 

%% 2.4: MANUAL Local streak removal and interpolation 
preprocessing_params.manualStreakRemoval_slices = 1:size(data_carried,3);
factor = 1;
preprocessing_params.manualStreakRemoval_factor = factor;
Y_local_removed = zeros(size(data_carried));

for s = preprocessing_params.manualStreakRemoval_slices
    [~, var_list, low_list] = streak_correction(data_carried(:,:,s),3,'plateau');
    figure; plot(low_list, var_list);
    [min_var,min_idx]=min(var_list);
    min_low = low_list(min_idx);

    [Y_local_removed(:,:,s), ~] = removeLocalStreaks(data_carried, s, factor*min_low, 3,'plateau');

    [Y_local_removed(:,:,s), ~] = interpolateLocalStreaks(Y_local_removed(:,:,s), 1, 0.8* min_low);
end
for s = 1:size(data_carried,3)
    if ~ismember(s, preprocessing_params.manualStreakRemoval_slices)
        Y_local_removed(:,:,s) = data_carried(:,:,s);
    end
end
data_carried = Y_local_removed;


%% 2.4: AUTO Local Streak Removal (No UI)
% This block runs streak correction, local streak removal, and interpolation in batch mode.

% define the slices to run 
preprocessing_params.autoStreakRemoval_slices = 50:100;
factor1 = 1;
factor2 = 1;
preprocessing_params.autoStreakRemoval_factor1 = factor1;
preprocessing_params.autoStreakRemoval_factor2 = factor2;
data_local_removed_auto = zeros(size(data_carried));
for s = preprocessing_params.autoStreakRemoval_slices
    % 1. Find optimal threshold automatically (e.g., by minimizing variance)
    [~, var_list, low_list] = streak_correction(data_carried(:,:,s), 3, 'valley');
    [~, min_idx] = min(var_list);
    min_low = low_list(min_idx);

    % 2. Remove local streaks automatically (no UI)
    [data_local_removed_auto(:,:,s), ~] = removeLocalStreaks(data_carried, s, factor1*min_low, 3, 'valley', true);

    % 3. Interpolate local streaks automatically (no UI)
    [data_local_removed_auto(:,:,s), ~] = interpolateLocalStreaks(data_local_removed_auto(:,:,s), 1, factor2*min_low, [], true);
    
    % print finished status
    fprintf('Finished slice %d\n', s);
end

% Copy back data_carried for slices not in autoStreakRemoval_slices
for s = 1:size(data_carried,3)
    if ~ismember(s, preprocessing_params.autoStreakRemoval_slices)
        data_local_removed_auto(:,:,s) = data_carried(:,:,s);
    end
end

data_carried = data_local_removed_auto;

%% 2.5 defect masking
% Normalize background 
[data_carried] = normalizeBackgroundToZeroMean3D(data_carried, rangetype, preprocessing_params.slice_normalize);

f1=figure;
d3gridDisplay(data_carried,'dynamic');
preprocessing_params.defect_slice = input('Enter defect slice number: ');
preprocessing_params.num_defect_type = input('enter how many types of defects to mask: ');
close(f1);
% methods: 
% 1. Gaussian window "gw"
% 2. truncated gaussian gaussian smoothing "tg"
% 3. thresholding and remove defect features "threshold"
preprocessing_params.defect_masking_method = 'tg';

switch preprocessing_params.defect_masking_method
    case 'gw'
        % Apply Gaussian window masking
        [data_masked, ~] = defect_masking(data_carried, preprocessing_params.defect_slice);
    case 'tg'
        % Apply flat disk mask with Gaussian smoothing
        % Interactive mask creation and application:
        %if isfield(preprocessing_params, 'defect_mask') && ~isempty(preprocessing_params.defect_mask)
            %[data_masked, ~] = gaussianMaskDefects(data_carried, [], [], preprocessing_params.defect_mask);
        %else
            [data_masked, preprocessing_params.defect_mask, defect_centers, sigmas] = gaussianMaskDefects(data_carried, preprocessing_params.defect_slice, preprocessing_params.num_defect_type);
        %end
    case 'threshold'
        % Apply threshold-based defect masking
        [data_masked, defect_mask] = thresholdDefects(data_carried, preprocessing_params.defect_slice);
    otherwise
        error('Unknown defect masking method. Choose "gw", "disk", or "threshold".');
end
data_carried = data_masked;

%% 2.6a: Correct streak
% Normalize background 
[data_carried] = normalizeBackgroundToZeroMean3D(data_carried, rangetype, preprocessing_params.slice_normalize);

[data_streakremoved, QPI_nostreaks] = RemoveStreaks(data_carried, 'Direction', 'vertical');
data_carried = data_streakremoved;

%% 2.6b: heal
% Normalize background 
[data_carried] = normalizeBackgroundToZeroMean3D(data_carried, rangetype, preprocessing_params.slice_normalize);

preprocessing_params.heal_direction = input('Enter direction to heal (horizontal/vertical/none): ', 's');

data_streakremoved_healed = heal_streaks(data_carried, preprocessing_params.heal_direction);

%% 2.6c: directional_plane (optional, zero slope at one direction)
% Normalize background 
[data_carried] = normalizeBackgroundToZeroMean3D(data_carried, rangetype, preprocessing_params.slice_normalize);

preprocessing_params.real_space_direction = 'horizontal';
[data_plane, mask] = d3plane_directional(data_carried, preprocessing_params.real_space_direction, 'LineWidth', 5);
data_carried = data_plane;

%% 2.end: Normalize background 
[Y] = normalizeBackgroundToZeroMean3D(data_carried, rangetype, preprocessing_params.slice_normalize);

%% 3. ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Save the preprocessed data~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
save('ZrSiTe_preprocessed0207_slice_1to80.mat', 'data_original', 'Y', 'data_cropped','data_masked',"energy_range", 'preprocessing_params')


%% Before Run Standardize
rangetype ='dynamic';
%% 4_pre Pick reference slice
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
window_type = {'gaussian', 5};
%window_type = '';

if same_size
    %[square_size] = squareDrawSize(Y_ref);
    square_size=[80,80];
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
%%
figure;
for k = 1:num_kernels
    subplot(1,num_kernels,k);
    imagesc(A1_ref{k});
    colorbar;
    axis square;
end

%% (ESS)noise level determination 
eta_data = estimate_noise(Y_ref, 'std');  

%% (Opt) determine SNR
SNR_data= var(A1{1,1}(:))/eta_data;
fprintf('SNR_data = %d', SNR_data);

%% Block 4: Find Optimal Activation for Reference Slice
% Set up display functions
figure;
dispfun = cell(1,  num_kernels);
for n = 1:num_kernels
    dispfun{n} = @(Y, A, X, kernel_sizes, kplus) showims(Y_ref, A1_ref{n}, X, A, X, kernel_sizes, kplus, 1);
end

% SBD settings.
miniloop_iteration = 1;
outerloop_maxIT= 8;
%params_ref.energy = energy_selected(params.ref_slice);
params_ref.lambda1 = [0.022, 0.020, 0.018, 0.02, 0.02];  % regularization parameter for Phase I
%params_ref.lambda1 = [0.15, 0.15, 0.15, 0.15, 0.15];  % regularization parameter for Phase I
params_ref.phase2 = false;
params_ref.kplus = ceil(0.2 * kernel_sizes);
params_ref.lambda2 = [0.03, 0.05, 0.05, 0.05];  % FINAL reg. param. value for Phase II
params_ref.nrefine = 4;
params_ref.signflip = 0.2;
params_ref.xpos = true;
params_ref.getbias = true;
params_ref.Xsolve = 'FISTA';

% noise variance for computeResidualQuality.m
params_ref.noise_var = eta_data;

% Update params for MTSBD
for k = 1:num_kernels
    params_ref.xinit{k}.X = X_ref(:,:,k);
    params_ref.xinit{k}.b = 0;
end

% Run and save
% 2. The fun part
[A_ref, X_ref, b_ref, extras_ref] = MT_SBD(Y_ref, kernel_sizes, params_ref, dispfun, A1_ref, miniloop_iteration, outerloop_maxIT);

%% Visualize Reference result 
[Y_rec,Y_rec_all] = visualizeRealResult(Y_ref,A_ref, X_ref, b_ref, extras_ref);


%% compare the Xout vs X manual
X0=zeros([size(Y_ref,1),size(Y_ref,2),length(defect_loc)]);
for i =1:length(defect_loc)
    X0(:,:,i)=locationsToMask(defect_loc{i},[size(Y_ref,1),size(Y_ref,2)]);
end 

[X_ref_aligned, ~, ~] = alignActivationMaps(X0, X_ref, kernel_sizes);
[X_similarity, ~] = computeActivationSimilarity(X0, X_ref_aligned, kernel_sizes,1);
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

params_ref.lambda1 = [0.03, 0.03, 0.03, 0.03];  % regularization parameter for Phase I
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
%% Reconstructed Y 
Y_rec_pad = zeros([size(Y_ref),num_kernels]);
for k = 1:num_kernels
    Y_rec_pad(:,:,k) = convfft2(A_ref_pad{1,k}, X_ref_pad(:,:,k)) + b_ref_pad(k);
end
%% Save the padded ones 
padfilename = sprintf('MTSBD_LiFeAs_%f meV.mat',1000*params_ref.energy);
%save(padfilename,'Y_ref','A_ref_pad', 'X_ref_pad', 'b_ref_pad', 'extras_ref_pad', 'params_ref');
save(padfilename,'Y_ref','A_ref', 'X_ref', 'b_ref', 'extras_ref', 'params_ref');

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
    axis square;
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
outerloop_maxIT= 3;

params.lambda1 = [0.179, 0.175, 0.171, 0.141, 0.141];  % regularization parameter for Phase I
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
    b_temp = 0; 
    params.xinit{k}.b = repmat(b_temp,[num_slices,1]);
end

if params.use_Xregulated
    [REG_Aout_ALL, REG_Xout_ALL, REG_bout_ALL, REG_extras_ALL] = MTSBD_Xregulated_all_slices(...
        Y_used, kernel_sizes, params, dispfun, A1_used, miniloop_iteration, outerloop_maxIT);
else
    [Aout_ALL, Xout_ALL, bout_ALL, ALL_extras] = MTSBD_all_slice(...
        Y_used, kernel_sizes, params, dispfun, A1_used, miniloop_iteration, outerloop_maxIT);
end

eta3dall = permute(repmat(eta_data3d,[outerloop_maxIT,1]),[2,1]);
observation_fidelity = eta3dall./squeeze(var(ALL_extras.phase1.residuals, 0, [1,2]));
save('ZrSiTe0207_meV_[80,80]_slice_1to80_new.mat', 'Y_used','Aout_ALL', 'Xout_ALL', 'bout_ALL', 'ALL_extras','params', 'eta_data3d','observation_fidelity');

% plot the observation fidelity  x axis is the number of slices, y is
% observation fidelity
figure; 
hold on 
for i = 1:outerloop_maxIT
    plot(1:num_slices,observation_fidelity(:,i));
    legend();
end
hold off
%% convert Aout_ALL to cell format
[num_slices, num_kernels] = size(bout_ALL);
Aout_ALL_cell = cell(num_slices, num_kernels);
for s = 1:num_slices
    for k = 1:num_kernels
        Aout_ALL_cell{s,k} = Aout_ALL{k}(:,:,s);
    end
end

%% Visualize Reference result 
for i = 40: 41
    pp=struct();
    pp.phase1.residuals = ALL_extras.phase1.residuals(:,:,i,:);
    pp.phase1.quality_metrics = ALL_extras.phase1.quality_metrics;
    visualizeRealResult(Y_used(:,:,i), Aout_ALL_cell(i,:), Xout_ALL, bout_ALL(i,:), pp);
end 
%% Kernels movies 
for k = 1:length(Aout_ALL)
    figure;
    d3gridDisplay(((Aout_ALL{k})), 'dynamic')
end
%% QPI movies 
for k = 1:length(Aout_ALL)
    figure;
    d3gridDisplay((qpiCalculate(Aout_ALL{k})), 'dynamic',-1)
end
%% Plot gaussianed activation
Xout_gaussian = zeros(size(Xout_ALL));
kernel_size = size(Aout_ALL_cell{1,1});
for i = 1:size(Xout_ALL, 3)
    Xout_gaussian(:,:,i) = Xout_gaussian_broaden(Xout_ALL(:,:,i), kernel_size);
    figure; imagesc(Xout_gaussian(:,:,i)); colormap("gray"); axis square
    title(sprintf('Kernel type %d', i))
end
%% Create reconstruction for all slices
Y_rec = zeros(size(Y_used));
for i = 1:size(Y_used,3)
    for k = 1:num_kernels
        Y_rec(:,:,i) = Y_rec(:,:,i) + convfft2(Aout_ALL_cell{i,k}, Xout_ALL(:,:,k)) + bout_ALL(i,k);
    end
end
figure;
d3gridDisplay(Y_rec, 'dynamic')

%% Create reconstruction for initialized all slices
Y_init = zeros(size(Y_used));
for i = 1:size(Y_used,3)
    for k = 1:num_kernels
        Y_init(:,:,i) = Y_init(:,:,i) + convfft2(A1_all{i,k}, X_ref(:,:,k));
    end
end
figure;
d3gridDisplay(Y_init, 'dynamic')

%% Create reconstruction for each kernel type
Y_rec_each = zeros([num_kernels,size(Y_used)]);
for i = 1:size(Y_used,3)
    for k = 1:num_kernels
        %Y_rec_each(k,:,:,i) = convfft2(Aout_ALL_cell{i,k}, Xout_ALL(:,:,k)) + bout_ALL(i,k);
        Y_rec_each(k,:,:,i) = convfft2(Aout_ALL_cell{i,k}, Xout_ALL(:,:,k));
    end
end

% create fft of Y_rec_each
FT_QPI_Y_rec_each = zeros([num_kernels,size(Y_used)]);
for k = 1:num_kernels
    FT_QPI_Y_rec_each(k,:,:,:) = qpiCalculate(squeeze(Y_rec_each(k,:,:,:)));
end

%% Normalize and combine Y_rec_each and its FT-QPI using method 2
% Reshape Y_rec_each and FT_QPI_Y_rec_each to combine all kernels
Y_rec_show_Full = [];
qpi_Y_rec_show_Full = [];
for k = 1:num_kernels
    Y_rec_show_Full = [Y_rec_show_Full, squeeze(Y_rec_each(k,:,:,:))];
    qpi_Y_rec_show_Full = [qpi_Y_rec_show_Full, squeeze(FT_QPI_Y_rec_each(k,:,:,:))];
end

% Normalize each slice across all kernels
for i = 1:size(Y_rec_show_Full,3)
    Y_rec_show_Full(:,:,i) = mat2gray(Y_rec_show_Full(:,:,i));
    qpi_Y_rec_show_Full(:,:,i) = 1-mat2gray(qpi_Y_rec_show_Full(:,:,i),[0,1]);
end

% Combine normalized reconstructions and their FT-QPI
Y_rec_ALL_show_norm = [Y_rec_show_Full; qpi_Y_rec_show_Full];

%% Display the normalized and combined results
figure;
d3gridDisplay(Y_rec_ALL_show_norm, 'dynamic');
title('Normalized Y_{rec}_-{each} and FT-QPI combined');

%% Y, Y_rec and Y residual

% Y and Y residual
Y_resi = ALL_extras.phase1.residuals(:,:,:,end);
Y_full_visualize = [Y_used, Y_rec, Y_resi];
qpi_Y_full_visualize = [qpiCalculate(Y_used), qpiCalculate(Y_rec), qpiCalculate(Y_resi)];

% Normalize each slice across all kernels
for i = 1:size(Y_full_visualize,3)
    Y_full_visualize(:,:,i) = mat2gray(Y_full_visualize(:,:,i));
    qpi_Y_full_visualize(:,:,i) = 1-mat2gray(qpi_Y_full_visualize(:,:,i),[0,1]);
end

% Combine normalized reconstructions and their FT-QPI
Y_show_norm = [Y_full_visualize; qpi_Y_full_visualize];

%% Display the Y_show_norm
figure;
d3gridDisplay(Y_show_norm, 'dynamic');
title('ALL combined');

%% 
Aout_show = [];
for i = 1: size(Aout_ALL,1)
    Aout_show = [Aout_show,mat2gray(Aout_ALL{i})];
end

figure;
d3gridDisplay(Aout_show, 'dynamic')

%% 
qpi_show = [];
for i = 1: size(Aout_ALL,1)
    qpi_show = [qpi_show,mat2gray(qpiCalculate(Aout_ALL{i}),[0,1])];
end

figure;
d3gridDisplay(qpi_show, 'dynamic')

%% Merge 3 ZrSiTe runs 
Aout_Full_energy = cell(2,1);
Aout_Full_energy{1,1}=C;
Aout_Full_energy{2,1}=D;

%%
save('ZrSiTe_kernel1&2_FULL_[80,80].mat', 'Y_used','Aout_Full_energy', 'Xout_A1', 'Xout_A2');

%% Show FULL
Aout_show_Full = [];
for i = 1: size(Aout_Full_energy,1)
    Aout_show_Full = [Aout_show_Full,Aout_Full_energy{i}];
end

figure;
d3gridDisplay(Aout_show_Full, 'dynamic')

%%
qpi_show_Full = [];
for i = 1: size(Aout_Full_energy,1)
    qpi_show_Full = [qpi_show_Full,qpiCalculate(Aout_Full_energy{i})];
end
figure;
d3gridDisplay(qpi_show_Full, 'dynamic',-1)

%%
Aout_show_norm = Aout_show_Full;
qpi_show_norm = qpi_show_Full;
for i = 1: 200
    Aout_show_norm(:,:,i) = mat2gray(Aout_show_Full(:,:,i));
    qpi_show_norm(:,:,i) = abs(mat2gray(qpi_show_Full(:,:,i),[0,1])-1);
end
ALL_show_norm = [Aout_show_norm;qpi_show_norm];
figure; 
d3gridDisplay(ALL_show_norm, 'dynamic');

%% write the video 
gridVideoWriter(rot90(Y_show_norm), V, 'dynamic', 100, 'invgray', 0, [800 800]);

% delay unit? 
%% QPI movies 
for k = 1:length(Aout_ALL)
    figure;
    d3gridDisplay(qpiCalculate(Aout_ALL{k}), 'dynamic')
end
%% Pad the output kernels to target sizes
target_size = [110, 110];  % Same target size as used for reference kernels
kernel_sizes_pad = repmat(target_size,[num_kernels,1]);

% Convert Aout_ALL to cell format with padding
A1_all_pad = cell(num_slices, num_kernels);
for s = 1:num_slices
    for k = 1:num_kernels
        % Get current kernel
        current_kernel = Aout_ALL_cell{s,k};
        sz = size(current_kernel);
        
        % Calculate padding needed
        pad_h = target_size(1) - sz(1);
        pad_w = target_size(2) - sz(2);
        
        % Calculate pre- and post-padding for centering
        pre_h = floor(pad_h / 2);
        post_h = ceil(pad_h / 2);
        pre_w = floor(pad_w / 2);
        post_w = ceil(pad_w / 2);
        
        % Pad the kernel to center it
        padded_kernel = padarray(current_kernel, [pre_h, pre_w], 'pre');
        padded_kernel = padarray(padded_kernel, [post_h, post_w], 'post');
        
        % Apply oblique projection
        A1_all_pad{s,k} = proj2oblique(padded_kernel);
    end
end

%% Convert A1_all_pad to matrix form and prepare noise
A1_all_pad_matrix = cell(num_kernels,1);
for k = 1:num_kernels
    A1_all_pad_matrix{k} = zeros(size(A1_all_pad{1,k},1),size(A1_all_pad{1,k},2),num_slices);
    for s = 1:num_slices
        A1_all_pad_matrix{k}(:,:,s) = A1_all_pad{s,k};
    end
end

eta_data3d = estimate_noise3D(Y, 'std');  

%% Run for all_pad_slice
% SBD settings.
miniloop_iteration = 5;
outerloop_maxIT= 5;

params.lambda1 = [5e-2, 5e-2,0.05];  % regularization parameter for Phase I
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

kernel_sizes_single = squeeze(max(kernel_sizes_pad,[],1));
Y_used = Y;
A1_used = A1_all_pad_matrix;

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
        Y_used, kernel_sizes_pad, params, dispfun, A1_used, miniloop_iteration, outerloop_maxIT);
else
    [Aout_ALL, Xout_ALL, bout_ALL, ALL_extras] = MTSBD_all_slice(...
        Y_used, kernel_sizes_pad, params, dispfun, A1_used, miniloop_iteration, outerloop_maxIT);
end

%% Save results
save('LiFeAs_slices.mat', 'Y_used', 'Aout_ALL', 'Xout_ALL', 'bout_ALL', 'ALL_extras', 'energy_selected');

%% convert Aout_ALL to cell format
Aout_ALL_cell = cell(num_slices, num_kernels);
for s = 1:num_slices
    for k = 1:num_kernels
        Aout_ALL_cell{s,k} = Aout_ALL{k}(:,:,s);
    end
end

%% Visualize result 
for i = 1: size(A1_all,1)
    pp=struct();
    pp.phase1.residuals = ALL_extras.phase1.residuals(:,:,i,:);
    pp.phase1.quality_metrics = ALL_extras.phase1.quality_metrics;
    visualizeRealResult(Y_used(:,:,i), Aout_ALL_cell(i,:), Xout_ALL, bout_ALL(i,:), pp);
end

%% Block 5: Sequential Processing with MT_SBD.m
% This block processes each slice individually using MT_SBD.m sequentially
fprintf('\n=== SEQUENTIAL PROCESSING BLOCK ===\n');
fprintf('Processing each slice individually with MT_SBD.m\n');

% Initialize storage for sequential results
Aout_seq = cell(num_slices, num_kernels);
Xout_seq = cell(num_slices, num_kernels);
bout_seq = zeros(num_slices, num_kernels);
extras_seq = cell(num_slices, 1);

% Estimate noise for all slices using estimate_noise3D (ROI defined once)
eta_data3d_seq = estimate_noise3D(Y, 'std');
fprintf('Noise estimation completed for all slices. ROI defined once.\n');

% Set up display functions (same as Block 4)
figure;
dispfun_seq = cell(1, num_kernels);

% SBD settings (same as Block 4)
miniloop_iteration_seq = 1;
outerloop_maxIT_seq = 8;

% Process each slice sequentially
for s = 1:num_slices
    fprintf('\n--- Processing Slice %d/%d ---\n', s, num_slices);
    slice_starttime = tic;

    for n = 1:num_kernels
        dispfun_seq{n} = @(Y, A, X, kernel_sizes, kplus) showims(Y, A1_all{s,n}, X, A, X, kernel_sizes, kplus, 1);
    end
    % Extract current slice
    Y_slice = Y(:,:,s);
    
    % Initialize kernels for this slice (same as Block 4)
    A1_slice = cell(1, num_kernels);
    for k = 1:num_kernels
        A1_slice{k} = A1_all{s,k};
    end
    
    % Set up parameters for this slice (same structure as Block 4)
    params_seq = struct();
    params_seq.lambda1 = [0.1, 0.2, 0.03, 0.03];  % Same as Block 4
    params_seq.phase2 = false;
    params_seq.kplus = ceil(0.2 * kernel_sizes);  % Same as Block 4
    params_seq.lambda2 = [0.05, 0.05, 0.05, 0.15];  % Same as Block 4
    params_seq.nrefine = 4;
    params_seq.signflip = 0.2;
    params_seq.xpos = true;
    params_seq.getbias = true;
    params_seq.Xsolve = 'FISTA';
    params_seq.noise_var = eta_data3d_seq(s);
    
    % Initialize xinit for this slice using reference results (same as Block 4)
    params_seq.xinit = cell(1, num_kernels);
    for k = 1:num_kernels
        params_seq.xinit{k}.X = X_ref(:,:,k);
        params_seq.xinit{k}.b = extras_ref.phase1.biter(k);
    end
    
    % Run MT_SBD for this slice (same as Block 4)
    try
        [Aout_slice, Xout_slice, bout_slice, extras_slice] = MT_SBD(...
            Y_slice, kernel_sizes, params_seq, dispfun_seq, A1_slice, ...
            miniloop_iteration_seq, outerloop_maxIT_seq);
        
        % Store results for this slice
        for k = 1:num_kernels
            Aout_seq{s,k} = Aout_slice{k};
            Xout_seq{s,k} = Xout_slice(:,:,k);
        end
        bout_seq(s,:) = bout_slice;
        extras_seq{s} = extras_slice;
        
        slice_runtime = toc(slice_starttime);
        fprintf('Slice %d completed in %.2fs\n', s, slice_runtime);
        
        % Save individual slice results immediately
        energy_val = energy_selected(s);
        % Format energy value for filename (convert to meV and handle negative values)
        energy_meV = energy_val * 1000;  % Convert to meV
        if energy_meV >= 0
            slice_filename = sprintf('LiFeAs_sequential_%+.1f_meV.mat', energy_meV);
        else
            slice_filename = sprintf('LiFeAs_sequential_%.1f_meV.mat', energy_meV);
        end
        slice_data = struct();
        slice_data.Aout = Aout_slice;
        slice_data.Xout = Xout_slice;
        slice_data.bout = bout_slice;
        slice_data.extras = extras_slice;
        slice_data.Y_slice = Y_slice;
        slice_data.slice_number = s;
        slice_data.energy = energy_val;
        slice_data.runtime = slice_runtime;
        slice_data.params = params_seq;
        
        save(slice_filename, '-struct', 'slice_data');
        fprintf('Saved individual results for slice %d (%.1f meV) to %s\n', s, energy_meV, slice_filename);
        
        % Optional: Visualize results for this slice
        if mod(s, 5) == 0 || s == num_slices  % Visualize every 5th slice and the last one
            fprintf('Visualizing results for slice %d...\n', s);
            visualizeRealResult(Y_slice, Aout_slice, Xout_slice, bout_slice, extras_slice);
        end
        
    catch ME
        fprintf('Error processing slice %d: %s\n', s, ME.message);
        % Save partial results even if there was an error
        if exist('Aout_slice', 'var') && exist('Xout_slice', 'var')
            energy_val = energy_selected(s);
            energy_meV = energy_val * 1000;  % Convert to meV
            if energy_meV >= 0
                slice_filename = sprintf('LiFeAs_sequential_%+.1f_meV_ERROR.mat', energy_meV);
            else
                slice_filename = sprintf('LiFeAs_sequential_%.1f_meV_ERROR.mat', energy_meV);
            end
            slice_data = struct();
            slice_data.Aout = Aout_slice;
            slice_data.Xout = Xout_slice;
            slice_data.bout = bout_slice;
            slice_data.extras = extras_slice;
            slice_data.Y_slice = Y_slice;
            slice_data.slice_number = s;
            slice_data.energy = energy_val;
            slice_data.error_message = ME.message;
            slice_data.params = params_seq;
            
            save(slice_filename, '-struct', 'slice_data');
            fprintf('Saved partial results for slice %d (%.1f meV, with error) to %s\n', s, energy_meV, slice_filename);
        end
        % Continue with next slice
        continue;
    end
end

% Save sequential results
save('ZrSiTe_sequential.mat', 'Aout_seq', 'Xout_seq', 'bout_seq', 'extras_seq', 'energy_selected');

% Convert sequential results to matrix format for comparison
Aout_seq_matrix = cell(num_kernels, 1);
for k = 1:num_kernels
    Aout_seq_matrix{k} = zeros(size(Aout_seq{1,k},1), size(Aout_seq{1,k},2), num_slices);
    for s = 1:num_slices
        if ~isempty(Aout_seq{s,k})
            Aout_seq_matrix{k}(:,:,s) = Aout_seq{s,k};
        end
    end
end

% Create kernel movies for sequential results
fprintf('\nCreating kernel movies for sequential results...\n');
for k = 1:num_kernels
    figure;
    d3gridDisplay(Aout_seq_matrix{k}, 'dynamic');
    title(sprintf('Sequential Kernel %d', k));
end

fprintf('\nSequential processing complete!\n');

%% Visualize sequential results
fprintf('\n=== VISUALIZING SEQUENTIAL RESULTS ===\n');

% Convert sequential results to matrix format for visualization
Aout_seq_matrix = cell(num_kernels, 1);
for k = 1:num_kernels
    Aout_seq_matrix{k} = zeros(size(Aout_seq{1,k},1), size(Aout_seq{1,k},2), num_slices);
    for s = 1:num_slices
        if ~isempty(Aout_seq{s,k})
            Aout_seq_matrix{k}(:,:,s) = Aout_seq{s,k};
        end
    end
end

% Calculate Y_reconstructed for each slice
Y_reconstructed = zeros(size(Y));
for s = 1:num_slices
    Y_slice_recon = zeros(size(Y,1), size(Y,2));
    for k = 1:num_kernels
        if ~isempty(Aout_seq{s,k}) && ~isempty(Xout_seq{s,k})
            Y_slice_recon = Y_slice_recon + convfft2(Aout_seq{s,k}, Xout_seq{s,k});
        end
    end
    Y_reconstructed(:,:,s) = Y_slice_recon;
end

% 1. Visualize optimized kernels (A_out) using d3gridDisplay
fprintf('Displaying optimized kernels (A_out) for each kernel type...\n');
for k = 1:num_kernels
    figure('Name', sprintf('Sequential Kernel %d (A_out)', k));
    d3gridDisplay(Aout_seq_matrix{k}, 'dynamic');
    title(sprintf('Sequential Processing - Kernel %d Evolution Across Slices', k));
    xlabel('X'); ylabel('Y'); zlabel('Slice/Energy');
end

% 2. Visualize reconstructed data (Y_reconstructed) using d3gridDisplay
fprintf('Displaying reconstructed data (Y_reconstructed)...\n');
figure('Name', 'Sequential Y_reconstructed');
d3gridDisplay(Y_reconstructed, 'dynamic');
title('Sequential Processing - Reconstructed Data Across Slices');
xlabel('X'); ylabel('Y'); zlabel('Slice/Energy');

% 3. Visualize original data for comparison
fprintf('Displaying original data for comparison...\n');
figure('Name', 'Original Data');
d3gridDisplay(Y, 'dynamic');
title('Original Data Across Slices');
xlabel('X'); ylabel('Y'); zlabel('Slice/Energy');

% 4. Visualize residual (difference between original and reconstructed)
fprintf('Displaying residual (Original - Reconstructed)...\n');
Y_residual = Y - Y_reconstructed;
figure('Name', 'Sequential Residual');
d3gridDisplay(Y_residual, 'dynamic');
title('Sequential Processing - Residual (Original - Reconstructed)');
xlabel('X'); ylabel('Y'); zlabel('Slice/Energy');

% 5. Calculate and display quality metrics
fprintf('Calculating quality metrics...\n');
residual_norm = zeros(num_slices, 1);
reconstruction_quality = zeros(num_slices, 1);

for s = 1:num_slices
    % Calculate residual norm for each slice
    residual_norm(s) = norm(Y_residual(:,:,s), 'fro') / norm(Y(:,:,s), 'fro');
    
    % Calculate reconstruction quality (1 - normalized residual)
    reconstruction_quality(s) = 1 - residual_norm(s);
end

% Plot quality metrics
figure('Name', 'Sequential Processing Quality Metrics');
subplot(2,1,1);
plot(energy_selected * 1000, residual_norm, 'b-o', 'LineWidth', 2);
title('Residual Norm vs Energy');
xlabel('Energy (meV)');
ylabel('Normalized Residual Norm');
grid on;

subplot(2,1,2);
plot(energy_selected * 1000, reconstruction_quality, 'r-o', 'LineWidth', 2);
title('Reconstruction Quality vs Energy');
xlabel('Energy (meV)');
ylabel('Reconstruction Quality');
grid on;

% Print summary statistics
fprintf('\n=== SEQUENTIAL PROCESSING QUALITY SUMMARY ===\n');
fprintf('Average Residual Norm: %.4f\n', mean(residual_norm));
fprintf('Average Reconstruction Quality: %.4f\n', mean(reconstruction_quality));
fprintf('Best Reconstruction Quality: %.4f (at %.1f meV)\n', max(reconstruction_quality), energy_selected(reconstruction_quality == max(reconstruction_quality)) * 1000);
fprintf('Worst Reconstruction Quality: %.4f (at %.1f meV)\n', min(reconstruction_quality), energy_selected(reconstruction_quality == min(reconstruction_quality)) * 1000);

fprintf('\nSequential results visualization complete!\n');

function Xout_gaussian = Xout_gaussian_broaden(X_out, kernel_sizes)
%XOUT_GAUSSIAN_BROADEN Apply Gaussian broadening to each channel of X_out based on kernel size.
%   Xout_gaussian = Xout_gaussian_broaden(X_out, kernel_sizes)
%   - X_out: [H x W x num_kernels] activation map
%   - kernel_sizes: [num_kernels x 2] array, each row is [height, width] for the kernel
%   - Xout_gaussian: [H x W x num_kernels] activation map after Gaussian broadening

    [h, w, num_kernels] = size(X_out);
    Xout_gaussian = zeros(h, w, num_kernels);
    for k = 1:num_kernels
        sigma = min(kernel_sizes(k,:)) / 10;
        window_size = ceil(3 * sigma);
        [x, y] = meshgrid(-window_size:window_size);
        gaussian_kernel = exp(-(x.^2 + y.^2)/(2*sigma^2));
        gaussian_kernel = gaussian_kernel / sum(gaussian_kernel(:));
        % Apply Gaussian broadening
        Xout_gaussian(:,:,k) = conv2(X_out(:,:,k), gaussian_kernel, 'same');
    end
end 