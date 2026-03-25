%% Centralized run configuration
cfg = init_config();

%% Before Run Standardize
rangetype ='dynamic';

%% Pick reference slice
figure;
d3gridDisplay(Y,rangetype);
params.ref_slice = input('Enter reference slice number: ');
num_kernels = input('enter the number of kernels you wish to apply: ');
if isempty(params.ref_slice) && ~isempty(cfg.reference.default_ref_slice)
    params.ref_slice = cfg.reference.default_ref_slice;
end
if isempty(num_kernels) && ~isempty(cfg.reference.default_num_kernels)
    num_kernels = cfg.reference.default_num_kernels;
end
close();

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
close();
%% Initialize reference kernels
% draw square on the data to include as many visible ripples of the scattering as possible 
same_size = 1;
kerneltype = "selected";
window_type = 'gaussian';
%window_type = '';

if same_size
    %[square_size] = squareDrawSize(Y_ref);
    square_size = [80,80];
    kernel_sizes = repmat(square_size,[num_kernels,1]);
    [A1_ref, A1_ref_crop] = initialize_kernels(Y_ref, num_kernels, kernel_sizes, kerneltype, window_type);
else
    A1_ref = cell(1, num_kernels);
    A1_ref_crop = cell(1, num_kernels);
    kernel_sizes = zeros(num_kernels, 2); % Store sizes of each kernel [height, width]
    for k = 1:num_kernels
        fprintf('Select region for kernel %d/%d\n', k, num_kernels);
        [square_size,position, mask] = squareDrawSize(Y_ref);           	% determine kernel size
        [A1_ref{k}, ~] = gridCropMask(Y_ref, mask);   % the cropped real data as kernel
        A1_ref_crop{k} = A1_ref{k};
        % Need to put each slice back onto the sphere
        A1_ref{k} = proj2oblique(A1_ref{k});
        % Store the kernel size
        kernel_sizes(k,:) = size(A1_ref_crop{k});
    end
end
for k = 1:num_kernels
    [A1_ref{k}, flipped_ref] = enforce_kernel_polarity(A1_ref{k}, A1_ref_crop{k});
    if flipped_ref
        fprintf('[kernel flip] reference kernel %d flipped\n', k);
    end
end
figure;
for k = 1:num_kernels
    subplot(1,num_kernels,k);
    imagesc(A1_ref{k});
    colorbar;
    axis square;
end

%% Noise level determination 
eta_data = estimate_noise(Y_ref, 'std');  

%% (Opt) determine SNR
SNR_data= var(A1_ref{1}(:))/eta_data;
fprintf('SNR_data = %d', SNR_data);

%% Find Optimal Activation for Reference Slice
% Set up display functions
figure;
dispfun = cell(1,  num_kernels);
for n = 1:num_kernels
    dispfun{n} = @(Y, A, X, kernel_sizes, kplus) showims(Y_ref, A1_ref{n}, X, A, X, kernel_sizes, kplus, 1);
end

% SBD settings.
miniloop_iteration = 1;
outerloop_maxIT= 2;
%params_ref.energy = energy_selected(params.ref_slice);
params_ref.lambda1 = [0.02, 0.02, 0.02, 0.02, 0.02];  % regularization parameter for Phase I
params_ref.phase2 = false;
params_ref.kplus = 0.2;
params_ref.lambda2 = [0.04, 0.04, 0.04, 0.04, 0.04];  % FINAL reg. param. value for Phase II
params_ref.nrefine = 4;
params_ref.signflip = 0.2;
params_ref.xpos = true;
params_ref.getbias = true;
params_ref.Xsolve = 'FISTA';

% noise variance for computeResidualQuality.m
params_ref.noise_var = eta_data;

% % Update params for MTSBD
% for k = 1:num_kernels
%     params_ref.xinit{k}.X = X_ref(:,:,k);
%     params_ref.xinit{k}.b = 0;
% end

% Run and save
% 2. The fun part
[A_ref, X_ref, b_ref, extras_ref] = MT_SBD(Y_ref, kernel_sizes, params_ref, dispfun, A1_ref, miniloop_iteration, outerloop_maxIT);

%% Visualize Reference result 
[Y_rec,Y_rec_all] = visualizeRealResult(Y_ref,A_ref, X_ref, b_ref, extras_ref);

%% Find Most Isolated Points and Initialize ALL Kernels (retire the most isolated points logic)
num_slices = size(Y,3);

% Choose method for kernel center selection
fprintf('Choose method for kernel center selection:\n');
fprintf('1. Find most isolated points automatically\n');
fprintf('2. Manually select 3 kernel centers\n');
choice = input('Enter choice (1 or 2): ');

if choice == 1
    fprintf('Calculating isolation scores and finding most isolated points...\n');
    
    % Choose target kernel sizes first
    type = cfg.isolation.target_kernel_size_type;
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
        threshold = max(X_ref(:,:,k),[],'all') / cfg.isolation.activation_threshold_divisor;
        %threshold = 0;
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
    title('Click on points to select kernel centers. Press Enter when done.');
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
A1_all_crop = cell(num_slices, num_kernels);
matrix = true;  % Use matrix form for A1
change_size = false;

for s = 1:num_slices
    fprintf('Initializing kernels for slice %d/%d...\n', s, num_slices);
    if matrix
        [A1_all(s,:), A1_all_crop(s,:)] = initialize_kernels_proliferation(Y(:,:,s), num_kernels, kernel_centers, window_type, target_kernel_sizes, 'interactive', change_size);
    else 
        [A1_all(s,:), A1_all_crop(s,:)] = initialize_kernels_proliferation(Y(:,:,s), num_kernels, kernel_centers, window_type, squeeze(kernel_sizes(s,:,:)), 'interactive', change_size);
    end
end
flip_slices_by_kernel = cell(1, num_kernels);
for s = 1:num_slices
    for k = 1:num_kernels
        [A1_all{s,k}, flipped] = enforce_kernel_polarity(A1_all{s,k}, A1_all_crop{s,k});
        if flipped
            flip_slices_by_kernel{k}(end+1) = s; %#ok<AGROW>
            fprintf('[kernel flip] kernel %d flipped at slice %d\n', k, s);
        end
    end
end
fprintf('==== Kernel flip summary (block init) ====\n');
for k = 1:num_kernels
    if isempty(flip_slices_by_kernel{k})
        fprintf('Kernel %d: flipped slices = []\n', k);
    else
        flip_slices_by_kernel{k} = unique(flip_slices_by_kernel{k});
        fprintf('Kernel %d: flipped slices = %s\n', k, mat2str(flip_slices_by_kernel{k}));
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

%% Orchastrate the truncation slices run 

% Configure the subset to run (global indices in the full dataset).
% Examples:
%   run_slice_idx = 1:80;
%   run_kernel_idx = [1,2,5];
run_slice_idx = 100:200;
run_kernel_idx = [1,2,3,4];

% Validate requested truncation indices.
if isempty(run_slice_idx) || any(run_slice_idx < 1) || any(run_slice_idx > size(Y,3))
    error('run_slice_idx must be non-empty and within 1:size(Y,3).');
end
if isempty(run_kernel_idx) || any(run_kernel_idx < 1) || any(run_kernel_idx > num_kernels)
    error('run_kernel_idx must be non-empty and within 1:num_kernels.');
end
run_slice_idx = unique(run_slice_idx(:).', 'stable');
run_kernel_idx = unique(run_kernel_idx(:).', 'stable');

% Build the truncated run inputs.
Y_used = Y(:,:,run_slice_idx);
A1_used = A1_all_matrix(run_kernel_idx);
for k = 1:numel(A1_used)
    A1_used{k} = A1_used{k}(:,:,run_slice_idx);
end
X_ref_used = X_ref(:,:,run_kernel_idx);
eta_data3d_used = eta_data3d(run_slice_idx);
kernel_sizes_used = kernel_sizes(run_kernel_idx, :);

% Keep run dimensions local to this orchestration.
run_num_slices = numel(run_slice_idx);
run_num_kernels = numel(run_kernel_idx);
fprintf('Truncated run configured: slices=%s, kernels=%s\n', ...
    mat2str(run_slice_idx), mat2str(run_kernel_idx));

%% initialize xinit for all slices with reference slice
params.xinit = cell(1, run_num_kernels);
for k = 1:run_num_kernels
    params.xinit{k}.X = X_ref_used(:,:,k);
    b_temp = 0; 
    params.xinit{k}.b = repmat(b_temp,[run_num_slices,1]);
end

%% Determine trusted-slice step weights (before all-slice run)
% Optional: trusted-slice weighting (can be disabled or left empty).
% Set this flag manually in this script if you want to use trusted slices
% (requires build_auto_trusted_slice_weights.m on the path).
use_trusted_weights = 0;           % 0 = unweighted (default), 1 = use trusted slices
trusted_ratio_threshold_default = 1.5;
use_default_manual_trusted_slices = true;
show_trusted_plot = true;

if use_trusted_weights
    % Automatically determine trusted slices for each kernel.
    % Criterion for slice i and kernel j:
    % std(A1_all_matrix{j}(:,:,i)) / std(noise_i) > threshold
    trusted_ratio_threshold = input(sprintf('Enter trusted-slice std-ratio threshold (e.g. %.2f): ', trusted_ratio_threshold_default));
    if isempty(trusted_ratio_threshold)
        trusted_ratio_threshold = trusted_ratio_threshold_default;
    end
    manual_trusted_slices = cell(1, run_num_kernels);
    if use_default_manual_trusted_slices && run_num_kernels >= 5
        manual_trusted_slices{1} = [1,4,5,8,10];
        manual_trusted_slices{2} = [1,5,8,9,10];
        manual_trusted_slices{3} = 7:11;
        manual_trusted_slices{4} = 7:11;
        manual_trusted_slices{5} = [3,5,6,7];
    end

    if exist('build_auto_trusted_slice_weights', 'file') ~= 2
        warning(['build_auto_trusted_slice_weights.m not found on path. ', ...
            'Falling back to unweighted mode (all slices treated as trusted).']);
        use_trusted_weights = 0;
    else
        [params.slice_weights, params.slice_weight_details] = ...
            build_auto_trusted_slice_weights(A1_used, eta_data3d_used, trusted_ratio_threshold, ...
            show_trusted_plot, manual_trusted_slices);
    end
end

if ~use_trusted_weights
    % Unweighted: allow MTSBD_all_slice to default slice_weights = ones
    params.slice_weights = [];
    params.slice_weight_details = struct();
    params.slice_weight_details.trusted_counts = run_num_slices * ones(1, run_num_kernels);
    params.slice_weight_details.method = "unweighted_all_trusted";
    params.slice_weight_details.trusted_ratio_threshold = trusted_ratio_threshold_default;
end

%% Run for all_slice
% SBD settings.
miniloop_iteration = 1;
outerloop_maxIT= 5;

lambda1_base_full = [0.020, 0.02, 0.02, 0.02, 0.018];
if numel(lambda1_base_full) < max(run_kernel_idx)
    error('lambda1_base_full must cover all selected kernel indices.');
end
params.lambda1_base = lambda1_base_full(run_kernel_idx);  % base regularization per selected kernel
trusted_counts = params.slice_weight_details.trusted_counts;
params.lambda1_weighted   = sqrt(trusted_counts) .* params.lambda1_base;
params.lambda1_unweighted = sqrt(run_num_slices) .* params.lambda1_base;
params.lambda1            = params.lambda1_unweighted;  % default/fallback for compatibility
%params.lambda1 = [0.15, 0.15, 0.15, 0.15, 0.15];  % regularization parameter for Phase I
params.phase2 = false;
params.kplus = ceil(0.2 * kernel_sizes_used);
lambda2_full = [0.04, 0.04, 0.04, 0.04, 0.04];
if numel(lambda2_full) < max(run_kernel_idx)
    error('lambda2_full must cover all selected kernel indices.');
end
params.lambda2 = lambda2_full(run_kernel_idx);  % FINAL reg. param. value for Phase II
params.nrefine = 4;
params.signflip = 0.2;
params.xpos = true;
params.getbias = true;
params.Xsolve =  'FISTA';
params.use_Xregulated = false;
params.noise_var = eta_data3d_used;
params.kernel_update_order = 1:run_num_kernels;  % local order for selected kernels

use_custom_order = true;
if use_custom_order
    use_custom_order = input('Use custom kernel update order for MTSBD_all_slice? (0/1): ');
end
if ~isempty(use_custom_order) && use_custom_order
    custom_order = input(sprintf('Enter kernel update permutation of 1:%d (e.g. [2 1 3 ...]): ', run_num_kernels));
    custom_order = custom_order(:).';
    if numel(custom_order) ~= run_num_kernels || any(custom_order < 1) || any(custom_order > run_num_kernels) || numel(unique(custom_order)) ~= run_num_kernels
        error('Invalid kernel update order. Must be a permutation of 1:run_num_kernels.');
    end
    params.kernel_update_order = custom_order;
end
fprintf('Kernel update order: %s\n', mat2str(params.kernel_update_order));

kernel_sizes_single = squeeze(max(kernel_sizes_used,[],1));
if use_trusted_weights
    fprintf('Trusted-slice weights ready. Counts per kernel: %s\n', mat2str(trusted_counts));
    fprintf('Lambda weighted (sqrt(trusted_count)): %s\n', mat2str(params.lambda1_weighted, 4));
else
    fprintf('Trusted-slice weighting disabled (unweighted mode; all slices treated as trusted).\n');
end
fprintf('Lambda unweighted (sqrt(total_slices)): %s\n', mat2str(params.lambda1_unweighted, 4));

% Set up display functions
figure;
dispfun = cell(1, run_num_kernels);
for n = 1:run_num_kernels
    dispfun{n} = @(Y, A, X, kernel_sizes_sing, kplus) showims(Y_used, A1_used{n}, X_ref_used(:,:,n), A, X, kernel_sizes_single, kplus, 1);
end

if params.use_Xregulated
    [REG_Aout_ALL, REG_Xout_ALL, REG_bout_ALL, REG_extras_ALL] = MTSBD_Xregulated_all_slices(...
        Y_used, kernel_sizes_used, params, dispfun, A1_used, miniloop_iteration, outerloop_maxIT);
else
    [Aout_ALL, Xout_ALL, bout_ALL, ALL_extras] = MTSBD_all_slice_modified(...
        Y_used, kernel_sizes_used, params, dispfun, A1_used, miniloop_iteration, outerloop_maxIT);
end

eta3dall = permute(repmat(eta_data3d_used,[outerloop_maxIT,1]),[2,1]);
observation_fidelity = eta3dall./squeeze(var(ALL_extras.phase1.residuals, 0, [1,2]));

% Save all-slice solver output. This legacy script does not use cfg, so we
% write to a local MAT file, avoiding overwriting an existing one.
if all(diff(run_slice_idx) == 1)
    slice_tag = sprintf('s%dto%d', run_slice_idx(1), run_slice_idx(end));
else
    slice_tag = sprintf('s%s', strrep(strrep(strrep(mat2str(run_slice_idx), ' ', ''), '[', ''), ']', ''));
end
kernel_tag = sprintf('k%s', strrep(strrep(strrep(mat2str(run_kernel_idx), ' ', ''), '[', ''), ']', ''));
allslice_file = sprintf('ZrSiTe0304_%s_%s_ALL.mat', slice_tag, kernel_tag);
if exist(allslice_file, 'file')
    ts = datestr(now, 'yyyymmdd_HHMMSS');
    [fpath, fname, fext] = fileparts(allslice_file);
    allslice_file = fullfile(fpath, sprintf('%s_%s%s', fname, ts, fext));
end
save(allslice_file, 'Y_used','Aout_ALL', 'Xout_ALL', 'bout_ALL', 'ALL_extras','A1_used','params', 'eta_data3d_used','observation_fidelity');
fprintf('Saved all-slice solver output to %s.\n', allslice_file);

% plot the observation fidelity  x axis is the number of slices, y is
% observation fidelity
figure; 
hold on 
for i = 1:outerloop_maxIT
    plot(1:run_num_slices,observation_fidelity(:,i));
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
%figure;
%d3gridDisplay(Y_rec, 'dynamic')

%% Create reconstruction for initialized all slices
Y_init = zeros(size(Y_used));

tic;
for i = 1:size(Y_used,3)
    for k = 1:num_kernels
        Y_init(:,:,i) = Y_init(:,:,i) + convfft2(A1_used{k}(:,:,i), X_ref_used(:,:,k));
    end
end
toc;

figure;
d3gridDisplay(Y_init, 'dynamic')

%% Reconstruction for per kernel initialize
Y_init_perkernel = zeros(size(Y_used,1),size(Y_used,2),size(Y_used,3),num_kernels);

tic;
for k = 1:num_kernels
    Y_init_perkernel(:,:,:,k) = convfft3(A1_used{k}, X_ref_used(:,:,k));
end

toc;
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
%Y_rec_ALL_show_norm = [Y_rec_show_Full; qpi_Y_rec_show_Full];

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

%% write the video 
gridVideoWriter(rot90(Y_show_norm), V, 'dynamic', 100, 'invgray', 0, [800 800]);

% delay unit? 
%% QPI movies 
for k = 1:length(Aout_ALL)
    figure;
    d3gridDisplay(qpiCalculate(Aout_ALL{k}), 'dynamic')
end

%% ~~~~~~~~~~~~~~~~~~~~~~~~~~Retired Blocks~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%% compare the Xout vs X manual (retire)
X0=zeros([size(Y_ref,1),size(Y_ref,2),length(defect_loc)]);
for i =1:length(defect_loc)
    X0(:,:,i)=locationsToMask(defect_loc{i},[size(Y_ref,1),size(Y_ref,2)]);
end 

[X_ref_aligned, ~, ~] = alignActivationMaps(X0, X_ref, kernel_sizes);
[X_similarity, ~] = computeActivationSimilarity(X0, X_ref_aligned, kernel_sizes,1);

%% Pad the A_ref to be size defined by user, normalize and use them as the A1 (retire)
target_size = cfg.sliceRunPadded.target_size;
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
%% Set ups before padded run (retire)
% Set up display functions
figure;
dispfun = cell(1, num_kernels);
for n = 1:num_kernels
    dispfun{n} = @(Y, A, X, kernel_sizes, kplus) showims(Y_ref, A1_ref{n}, X, A, X, kernel_sizes, kplus, 1);
end

% SBD settings.
miniloop_iteration = cfg.sliceRunPadded.miniloop_iteration;
outerloop_maxIT= cfg.sliceRunPadded.outerloop_maxIT;

params_ref.lambda1 = cfg.sliceRunPadded.lambda1;  % regularization parameter for Phase I
%params_ref.lambda1 = [0.15, 0.15, 0.15, 0.15, 0.15];  % regularization parameter for Phase I
params_ref.phase2 = cfg.sliceRunPadded.phase2;
params_ref.kplus = ceil(cfg.sliceRunPadded.kplus_factor * kernel_sizes);
params_ref.lambda2 = cfg.sliceRunPadded.lambda2;  % FINAL reg. param. value for Phase II
params_ref.nrefine = cfg.sliceRunPadded.nrefine;
params_ref.signflip = cfg.sliceRunPadded.signflip;
params_ref.xpos = cfg.sliceRunPadded.xpos;
params_ref.getbias = cfg.sliceRunPadded.getbias;
params_ref.Xsolve = cfg.sliceRunPadded.Xsolve;

% noise variance for computeResidualQuality.m
params_ref.noise_var = eta_data;
%% Run the padded initialization (retire)
% 2. The fun part
[A_ref_pad, X_ref_pad, b_ref_pad, extras_ref_pad] = MT_SBD(Y_ref, kernel_sizes_pad, params_ref, dispfun, A1_ref, miniloop_iteration, outerloop_maxIT);
%% Visualize Padded result (retire)
visualizeRealResult(Y_ref,A_ref_pad, X_ref_pad, b_ref_pad, extras_ref_pad);
%% Reconstructed Y (retire) 
Y_rec_pad = zeros([size(Y_ref),num_kernels]);
for k = 1:num_kernels
    Y_rec_pad(:,:,k) = convfft2(A_ref_pad{1,k}, X_ref_pad(:,:,k)) + b_ref_pad(k);
end
%% Save the padded ones (retire) 
padfilename = sprintf('MTSBD_LiFeAs_%f meV.mat',1000*params_ref.energy);
%save(padfilename,'Y_ref','A_ref_pad', 'X_ref_pad', 'b_ref_pad', 'extras_ref_pad', 'params_ref');
save(padfilename,'Y_ref','A_ref', 'X_ref', 'b_ref', 'extras_ref', 'params_ref');

%% Insert previous results (retire)

for i = 1:5
    % adjust the inplane shift
    A1_all_matrix{i} = inplaneShift(A1_all_matrix{i},[41,41],[40,41]);
    A1_all_matrix{i}(:,:,80:130)=Aout_ALL{i};
end

%% (retire)
for i = 1:5
    % adjust the inplane shift
    A1_all_matrix{i} = A1_all_matrix{i}(:,:,70:90);
end

%% Visualize Reference result (retire)
for i = 40: 41
    pp=struct();
    pp.phase1.residuals = ALL_extras.phase1.residuals(:,:,i,:);
    pp.phase1.quality_metrics = ALL_extras.phase1.quality_metrics;
    visualizeRealResult(Y_used(:,:,i), Aout_ALL_cell(i,:), Xout_ALL, bout_ALL(i,:), pp);
end 
%% Show FULL&QPI&QPI_padded (retire)
Aout_show_Full = [];
for i = 1: size(Aout_Full_energy,1)
    Aout_show_Full = [Aout_show_Full,Aout_Full_energy{i}];
end

figure;
d3gridDisplay(Aout_show_Full, 'dynamic')

qpi_show_Full = [];
for i = 1: size(Aout_Full_energy,1)
    qpi_show_Full = [qpi_show_Full,qpiCalculate(Aout_Full_energy{i})];
end
figure;
d3gridDisplay(qpi_show_Full, 'dynamic',-1)
%% (retire)
conv2
qpi_show_Full_pad = [];
for i = 1: size(Aout_Full_energy,1)
    mid=qpiCalculate(Aout_Full_energy{i},496);
    mid = (mid-min(mid,[],'all'))/max(mid,[],'all');
    qpi_show_Full_pad = [qpi_show_Full_pad,mid];
end
mid = qpiCalculate(Y_used);
qpi_show_Full_pad = [qpi_show_Full_pad, (mid-min(mid,[],'all'))/max(mid,[],'all')];
figure;
d3gridDisplay(qpi_show_Full_pad, 'dynamic',-1)
%% (retire)
Aout_show_norm = Aout_show_Full;
qpi_show_norm = qpi_show_Full;
for i = 1: 200
    Aout_show_norm(:,:,i) = mat2gray(Aout_show_Full(:,:,i));
    qpi_show_norm(:,:,i) = abs(mat2gray(qpi_show_Full(:,:,i),[0,1])-1);
end
ALL_show_norm = [Aout_show_norm;qpi_show_norm];
figure; 
d3gridDisplay(ALL_show_norm, 'dynamic');

%% Pad the output kernels to target sizes (retire)
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

%% Convert A1_all_pad to matrix form and prepare noise (retire)
A1_all_pad_matrix = cell(num_kernels,1);
for k = 1:num_kernels
    A1_all_pad_matrix{k} = zeros(size(A1_all_pad{1,k},1),size(A1_all_pad{1,k},2),num_slices);
    for s = 1:num_slices
        A1_all_pad_matrix{k}(:,:,s) = A1_all_pad{s,k};
    end
end

eta_data3d = estimate_noise3D(Y, 'std');  

%% Run for all_pad_slice (retire)
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

%% Save results (retire)
save('LiFeAs_slices.mat', 'Y_used', 'Aout_ALL', 'Xout_ALL', 'bout_ALL', 'ALL_extras', 'energy_selected');

%% convert Aout_ALL to cell format (retire)
Aout_ALL_cell = cell(num_slices, num_kernels);
for s = 1:num_slices
    for k = 1:num_kernels
        Aout_ALL_cell{s,k} = Aout_ALL{k}(:,:,s);
    end
end

%% Visualize result (retire) 
for i = 1: size(A1_all,1)
    pp=struct();
    pp.phase1.residuals = ALL_extras.phase1.residuals(:,:,i,:);
    pp.phase1.quality_metrics = ALL_extras.phase1.quality_metrics;
    visualizeRealResult(Y_used(:,:,i), Aout_ALL_cell(i,:), Xout_ALL, bout_ALL(i,:), pp);
end

%% Block 5: Sequential Processing with MT_SBD.m (retire)
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

%% Display the Y_show_norm (retire)
figure;
d3gridDisplay(Y_show_norm(), 'dynamic');
title('ALL combined');

%%  (retire)
Aout_show = [];
for i = 1: size(Aout_ALL,1)
    Aout_show = [Aout_show,mat2gray(Aout_ALL{i})];
end

figure;
d3gridDisplay(Aout_show, 'dynamic')

%%  (retire)
qpi_show = [];
for i = 1: size(Aout_ALL,1)
    qpi_show = [qpi_show,mat2gray(qpiCalculate(Aout_ALL{i}),[0,1])];
end

figure;
d3gridDisplay(qpi_show, 'dynamic')

%% Merge 3 ZrSiTe runs (retire) 
Aout_Full_energy = cell(2,1);
Aout_Full_energy{1,1}=C;
Aout_Full_energy{2,1}=D;
save('ZrSiTe_kernel1&2_FULL_[80,80].mat', 'Y_used','Aout_Full_energy', 'Xout_A1', 'Xout_A2');


%% Visualize sequential results (retire)
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

function [slice_weights, details] = build_auto_trusted_slice_weights(A1_all_matrix, noise_var, threshold, show_plot, manual_trusted_slices)
%BUILD_AUTO_TRUSTED_SLICE_WEIGHTS Build binary trusted-slice weights per kernel.
%   Inputs:
%       A1_all_matrix : cell array, each cell is [h x w x num_slices] kernel stack
%       noise_var     : scalar or [num_slices x 1] noise variance
%       threshold     : trusted ratio threshold
%       show_plot     : whether to show ratio-vs-slice visualization
%       manual_trusted_slices : cell array of manual trusted slice indices (optional)
%   Rule:
%       ratio(i,j) = std(A1_all_matrix{j}(:,:,i)) / sqrt(noise_var(i))
%       trusted if ratio(i,j) > threshold

    if nargin < 4 || isempty(show_plot)
        show_plot = true;
    end
    if nargin < 5 || isempty(manual_trusted_slices)
        manual_trusted_slices = {};
    end

    num_kernels = numel(A1_all_matrix);
    num_slices = size(A1_all_matrix{1}, 3);
    eps0 = 1e-12;

    if isscalar(noise_var)
        noise_var = repmat(noise_var, [num_slices, 1]);
    else
        noise_var = noise_var(:);
    end
    if numel(noise_var) ~= num_slices
        error('noise_var must be scalar or one value per slice.');
    end

    noise_std = sqrt(max(double(noise_var), eps0));
    ratio = zeros(num_slices, num_kernels);
    slice_weights = zeros(num_slices, num_kernels);
    trusted_slices = cell(1, num_kernels);
    trusted_counts = zeros(1, num_kernels);

    for k = 1:num_kernels
        Ak = A1_all_matrix{k};
        if size(Ak, 3) ~= num_slices
            error('All kernel stacks in A1_all_matrix must have the same num_slices.');
        end

        for s = 1:num_slices
            curr_slice = double(Ak(:,:,s));
            ratio(s,k) = std(curr_slice(:), 0) / noise_std(s);
        end

        idx = find(ratio(:,k) > threshold);
        if isempty(idx)
            % Ensure at least one trusted slice: keep the max-ratio slice.
            [~, idx_max] = max(ratio(:,k));
            idx = idx_max;
        end

        slice_weights(idx, k) = 1;
        trusted_slices{k} = idx(:).';
        trusted_counts(k) = numel(idx);
    end

    if show_plot
        figure('Name', 'Trusted Slice Ratio by Kernel');
        hold on;
        h = gobjects(num_kernels,1);
        for k = 1:num_kernels
            h(k) = scatter(1:num_slices, ratio(:,k), 18, 'filled');
            if numel(manual_trusted_slices) >= k && ~isempty(manual_trusted_slices{k})
                idxm = unique(round(manual_trusted_slices{k}(:).'));
                idxm = idxm(idxm >= 1 & idxm <= num_slices);
                if ~isempty(idxm)
                    c = h(k).CData;
                    scatter(idxm, ratio(idxm,k), 60, c, 'o', 'LineWidth', 1.2);
                end
            end
        end
        yline(threshold, 'k--', 'LineWidth', 1.5);
        xlabel('Slice index');
        ylabel('std(kernel slice) / std(noise)');
        title('Trusted-slice std ratio vs slice index');
        legend_labels = arrayfun(@(k) sprintf('Kernel %d', k), 1:num_kernels, 'UniformOutput', false);
        h_manual = scatter(nan, nan, 60, 'o', 'k', 'LineWidth', 1.2);
        h_thr = plot(nan, nan, 'k--', 'LineWidth', 1.5);
        legend([h; h_manual; h_thr], [legend_labels, {'Manual trusted slices', 'Threshold'}], 'Location', 'best');
        grid on;
        hold off;
    end

    details = struct();
    details.threshold = threshold;
    details.ratio = ratio;
    details.trusted_slices = trusted_slices;
    details.trusted_counts = trusted_counts;
end

function [A_out, flipped] = enforce_kernel_polarity(A_in, A_anchor)
    % Resolve sign ambiguity after kernel normalization/projection.
    % Without this, equivalent kernels can differ by a global +/-1 factor.
    A_out = A_in;
    flipped = false;
    if nargin < 2 || isempty(A_anchor)
        [~, idx] = max(abs(A_in(:)));
        score = A_in(idx);
    else
        score = sum(A_in(:) .* A_anchor(:));
    end
    if score < 0
        A_out = -A_in;
        flipped = true;
    end
end