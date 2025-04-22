%% SBD with Proximity Regularization - Script Version
% This script implements blind deconvolution with proximity regularization
% Each block can be run independently for testing and debugging

%% Block 1: Generate Test Data and Initialize Kernels
% Generate test set
SNR = 2;              % Signal-to-noise ratio
N_obs = 50;           % Observation lattice size (50x50)
observation_resolution = 4;  % Resolution: 3 pixels per lattice site
defect_density = 0.005;      % 0.1% defect density
num_slices = 10;
% Generate data with custom parameters
[Y, A0, X0, params] = properGen_full(SNR, N_obs, observation_resolution, defect_density, ...
    'LDoS_path', 'example_data/LDoS_single_defect_save.mat', ...
    'num_slices', num_slices);

num_kernels = params.num_kernels;
num_slices = params.num_slices; 
kernel_sizes_all = cell(num_slices,1);
for s = 1: num_slices
    kernel_sizes_all{s} = zeros(num_kernels, 2);
    for k = 1: num_kernels
    kernel_sizes_all{s}(k,:) = size(A0{s,k}); 
    end
end

% Normalize and visualize
rangetype = 'dynamic';  
Y = normalizeBackgroundToZeroMean3D(Y, rangetype); 
Y = proj2oblique(Y);

% Display the 3D data and let the user input the reference slice
figure; 
d3gridDisplay(Y,'dynamic');
% Let the user input the reference slice
title('Select a reference slice');
params.ref_slice = input('Enter reference slice number: ');

% Validate the input
if isempty(params.ref_slice) || ~isnumeric(params.ref_slice) || params.ref_slice < 1 || params.ref_slice > size(Y, 3)
    error('Invalid reference slice. Please enter a number between 1 and %d.', size(Y, 3));
end

fprintf('Using slice %d as reference\n', params.ref_slice);

Y_ref = Y(:,:,params.ref_slice);
X0_ref = X0;
A0_ref = cell(1, params.num_kernels);

for i= 1: params.num_kernels
    A0_ref{i} = A0{params.ref_slice,i};
end
% Display the selected reference slice
figure;
imagesc(Y_ref);
colorbar;
title(sprintf('Reference Slice %d', params.ref_slice));
axis square;
kernel_size_ref = squeeze(params.kernel_sizes(params.ref_slice,:,:));

% Initialize kernels
kerneltype = 'selected';
window_type = {'gaussian', 5};  % Example: gaussian window with alpha
A1 = initialize_kernels(Y_ref, params.num_kernels, kernel_size_ref, kerneltype, window_type);

% Display initialized kernels
figure;
for n = 1:params.num_kernels
    subplot(1, params.num_kernels, n);
    imagesc(A1{n});
    title(sprintf('Initial Kernel %d', n));
    colorbar;
    axis square;
end

%% Block 2: Find Optimal Activation for Reference Slice
% Set up display functions
figure;
dispfun = cell(1, params.num_kernels);
for n = 1:num_kernels
    dispfun{n} = @(Y, A, X, kernel_sizes, kplus) showims(Y_ref, A1{n}, X0_ref(:,:,n), A, X, kernel_sizes, kplus, 1);
end

% SBD settings
initial_iteration = 3;
maxIT = 30;
params.phase2 = false;
params.lambda1 = [5e-2, 5e-2, 5e-2];  % regularization parameter for Phase I
params.kplus = ceil(0.5 * params.kernel_sizes);
params.lambda2 = [1e-2, 1e-2, 1e-2];  % FINAL reg. param. value for Phase II
params.signflip = 0.2;
params.xpos = true;
params.getbias = true;
params.Xsolve = 'FISTA';
params.Xsolve = 'FISTA';
% Store ground truth for testing
params.X0 = X0_ref;
params.A0 = A0_ref;
%{
for k = 1:params.num_kernels
    params.xinit{k}.X = X0_ref(:,:,k);
    params.xinit{k}.b = extras.phase1.biter(k); 
end
%}

% Run SBD on reference slice
fprintf('Finding optimal activation for reference slice %d...\n', params.ref_slice);
[A_ref, X_ref, ~, extras] = SBD_test_multi_demixing(...
    Y_ref, kernel_size_ref, params, dispfun, A1, initial_iteration, maxIT);
% Store reference results
ref_results = struct();
ref_results.A = A_ref;
ref_results.X = X_ref;
ref_results.extras = extras;

%% Block 3: Find Most Isolated Points and Generate Initial A Guesses
% This block calculates isolation scores for defects and generates initial kernel guesses

fprintf('Calculating isolation scores and finding most isolated points...\n');

% Get defect positions from X_ref
most_isolated_points = cell(1, params.num_kernels);
isolation_scores = cell(1, params.num_kernels);

% First, get all defect positions for each kernel
defect_positions = cell(1, params.num_kernels);
num_defects = zeros(1, params.num_kernels);

% Create figure for activation value distributions
figure('Name', 'Activation Value Distributions');
for k = 1:params.num_kernels
    % Plot histogram of activation values
    subplot(2, params.num_kernels, k);
    activation_values = X_ref(:,:,k);
    h = histogram(activation_values(activation_values > 0), 50);
    set(gca, 'YScale', 'log');  % Set y-axis to log scale
    title(sprintf('Kernel %d Activation Distribution', k));
    xlabel('Activation Value');
    ylabel('Frequency (log scale)');
    
    % Add vertical line for threshold
    threshold = max(X_ref(:,:,k),[],'all')/10;
    hold on;
    xline(threshold, 'r--', 'Threshold');
    hold off;
    
    % Plot cumulative distribution
    subplot(2, params.num_kernels, k + params.num_kernels);
    [counts, edges] = histcounts(activation_values(activation_values > 0), 50, 'Normalization', 'cdf');
    stairs(edges(1:end-1), counts);
    title(sprintf('Kernel %d Cumulative Distribution', k));
    xlabel('Activation Value');
    ylabel('Cumulative Frequency');
    hold on;
    xline(threshold, 'r--', 'Threshold');
    hold off;
    
    % Get positions of defects above threshold
    [rows, cols] = find(X_ref(:,:,k) > threshold);  % Note: X_ref is [spatial_x, spatial_y, num_kernels]
    defect_positions{k} = [rows, cols];
    num_defects(k) = size(defect_positions{k}, 1);
    fprintf('Kernel %d has %d defects above threshold %.4f\n', k, num_defects(k), threshold);
end

% Calculate isolation scores for each kernel
for k = 1:params.num_kernels
    if num_defects(k) == 0
        warning('No defects found for kernel %d', k);
        continue;
    end
    
    % Create summed activation map of all other kernels
    X_others = zeros(size(X_ref(:,:,1)));  % Initialize with correct spatial dimensions
    for l = 1:params.num_kernels
        if l ~= k
            X_others = X_others + X_ref(:,:,l);  % Sum all other kernels
        end
    end
    
    % Get positions of defects in other kernels
    [other_rows, other_cols] = find(X_others > max(X_others,[],'all')/10);
    other_positions = [other_rows, other_cols];
    
    % Get default kernel size for boundary check
    default_kernel_size = size(A_ref{k});
    half_kernel_size = floor(default_kernel_size/2);
    
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
        % Calculate distances to all points in X_others
        diffs = other_positions - valid_defects(i,:);
        distances = sum(diffs.^2, 2);  % squared Euclidean distances
        S_k(i) = min(distances);  % minimum distance to any defect in other kernels
    end
    
    % Find the most isolated point among valid points
    [max_score, max_idx] = max(S_k);
    
    % Convert back to original indices
    valid_indices = find(valid_points);
    max_idx = valid_indices(max_idx);
    
    most_isolated_points{k} = defect_positions{k}(max_idx,:);
    isolation_scores{k} = S_k;
    
    fprintf('Kernel %d: Most isolated point at (%d,%d) with score %.2f\n', ...
        k, most_isolated_points{k}(1), most_isolated_points{k}(2), max_score);
end

% Visualize the results
figure('Name', 'Isolation Analysis');
for k = 1:params.num_kernels
    subplot(2, params.num_kernels, k);
    
    % Plot activation map for current kernel
    imagesc(X_ref(:,:,k));  % Directly index the k-th kernel
    hold on;
    
    % Plot all defects
    scatter(defect_positions{k}(:,2), defect_positions{k}(:,1), 50, 'w', 'o');
    
    % Highlight the most isolated point
    if ~isempty(most_isolated_points{k})
        scatter(most_isolated_points{k}(2), most_isolated_points{k}(1), 100, 'r', '*');
    end
    
    title(sprintf('Kernel %d', k));
    colorbar;
    axis square;
    hold off;
    
    % Plot summed map of other kernels
    subplot(2, params.num_kernels, k + params.num_kernels);
    X_others = zeros(size(X_ref(:,:,1)));
    for l = 1:params.num_kernels
        if l ~= k
            X_others = X_others + X_ref(:,:,l);
        end
    end
    imagesc(X_others);
    hold on;
    % Highlight the most isolated point from kernel k on the summed map
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

fprintf('Isolation analysis complete.\n');

%% Block 4: Initialize Kernels Using Most Isolated Points
use_ref = 1;
change_size = false;

% This block uses the most isolated points to generate initial kernel guesses for all slices

fprintf('\nInitializing kernels using most isolated points for all slices...\n');

% Initialize arrays for reference kernel sizes
ref_kernel_sizes = zeros(params.num_kernels, 2);
for k = 1:params.num_kernels
    ref_kernel_sizes(k,:) = size(A_ref{k});
end

% Use kernel_sizes_all directly as target kernel sizes
target_kernel_sizes = kernel_sizes_all;

% align the most isolated points with ground truth activation 
[~, offset, ~] = alignActivationMaps(X0_ref, X_ref, kernel_size_ref);
% apply this offset to all most isolated points
for k = 1:params.num_kernels
    most_isolated_points{k} = most_isolated_points{k} + offset(k,:);
end

% display the most isolated points on the ground truth activation map and let the user adjust the Points
figure;
imagesc(Y_ref);
hold on;
for k = 1:params.num_kernels
    scatter(most_isolated_points{k}(2), most_isolated_points{k}(1), 100, 'r', '*');
end

% Convert most_isolated_points from cell to matrix format
kernel_centers = zeros(params.num_kernels, 2);
for k = 1:params.num_kernels
    if ~isempty(most_isolated_points{k})
        kernel_centers(k,:) = most_isolated_points{k};
    else
        error('No isolated point found for kernel %d', k);
    end
end

% Initialize kernels for all slices
A1_all = cell(params.num_slices, params.num_kernels);

for s = 1:params.num_slices
    fprintf('Initializing kernels for slice %d/%d...\n', s, params.num_slices);
    if use_ref
    % Initialize kernels with ref sizes for this slice
    A1_all(s,:) = initialize_kernels_proliferation(Y(:,:,s), params.num_kernels, kernel_centers, window_type, ref_kernel_sizes, 'interactive', change_size);
    else 
    % Initialize kernels with target sizes for this slice
    A1_all(s,:) = initialize_kernels_proliferation(Y(:,:,s), params.num_kernels, kernel_centers, window_type, target_kernel_sizes{s}, 'interactive', change_size);
    end
end

% Visualize initialized kernels for reference slice
figure('Name', 'Initialized Kernels (Reference Slice)');
for k = 1:params.num_kernels
    subplot(1,params.num_kernels,k);
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
init_results.kernel_sizes = ref_kernel_sizes;

fprintf('Kernel initialization complete for all slices.\n');
%% Test block for perform SBD on all slices at once  

% Initialize arrays for results (these will be written to in parallel)
A_all_out = cell(params.num_slices, params.num_kernels);
X_all_out = zeros(size(Y,1), size(Y,2), params.num_kernels, params.num_slices);
b_all_out = zeros(params.num_kernels, params.num_slices);
extras_all_out = cell(params.num_slices, 1);


    Y_used = Y(:,:,slice);
    A1_used = A1_all(slice,:);
    kernel_sizes_used=zeros(size(A1_used,2),2);

    for i = 1:size(A1_used,2)
        kernel_sizes_used(i,:) = size(A1_used{i});
    end
    % Set up display functions
    figure;
    dispfun = cell(1, params.num_kernels);
    for n = 1:params.num_kernels
        dispfun{n} = @(Y, A, X, kernel_sizes, kplus) showims(Y_used, A1_used{n}, X0_ref(:,:,n), A, X, kernel_sizes, kplus, 1);
    end

    % SBD settings
    initial_iteration = 3;
    maxIT = 30;
    params.phase2 = false;
    params.lambda1 = [5e-2, 5e-2, 5e-2];  % regularization parameter for Phase I
    params.kplus = ceil(0.5 * params.kernel_sizes);
    params.lambda2 = [1e-2, 1e-2, 1e-2];  % FINAL reg. param. value for Phase II
    params.signflip = 0.2;
    params.xpos = true;
    params.getbias = true;
    params.Xsolve = 'FISTA';
    params.Xsolve = 'FISTA';
    % Store ground truth for testing
    params.X0 = X0_ref;
    params.A0 = A1_used;

    for k = 1:params.num_kernels
        params.xinit{k}.X = X_ref(:,:,k);
        params.xinit{k}.b = extras.phase1.biter(k); 
    end


    % Run SBD on reference slice
    fprintf('Finding optimal activation for reference slice %d...\n', params.ref_slice);
    [Aout_slice, Xout_slice, bout_slice, slice_extras] = SBD_test_multi_demixing(...
        Y_used, kernel_sizes_used, params, dispfun, A1_used, initial_iteration, maxIT);

    A_all_out(slice,:) = Aout_slice;
    X_all_out(:,:,:,slice) = Xout_slice;
    b_all_out(:,slice) = bout_slice;
    extras_all_out{slice} = slice_extras;

    fprintf('Completed slice %d/%d\n', slice, params.num_slices);

%% Test block for xinit, use slice 

% Initialize arrays for results (these will be written to in parallel)
A_all_out = cell(params.num_slices, params.num_kernels);
X_all_out = zeros(size(Y,1), size(Y,2), params.num_kernels, params.num_slices);
b_all_out = zeros(params.num_kernels, params.num_slices);
extras_all_out = cell(params.num_slices, 1);

for slice=1:params.num_slices
    Y_used = Y(:,:,slice);
    A1_used = A1_all(slice,:);
    kernel_sizes_used=zeros(size(A1_used,2),2);

    for i = 1:size(A1_used,2)
        kernel_sizes_used(i,:) = size(A1_used{i});
    end
    % Set up display functions
    figure;
    dispfun = cell(1, params.num_kernels);
    for n = 1:params.num_kernels
        dispfun{n} = @(Y, A, X, kernel_sizes, kplus) showims(Y_used, A1_used{n}, X0_ref(:,:,n), A, X, kernel_sizes, kplus, 1);
    end

    % SBD settings
    initial_iteration = 3;
    maxIT = 30;
    params.phase2 = false;
    params.lambda1 = [5e-2, 5e-2, 5e-2];  % regularization parameter for Phase I
    params.kplus = ceil(0.5 * params.kernel_sizes);
    params.lambda2 = [1e-2, 1e-2, 1e-2];  % FINAL reg. param. value for Phase II
    params.signflip = 0.2;
    params.xpos = true;
    params.getbias = true;
    params.Xsolve = 'FISTA';
    params.Xsolve = 'FISTA';
    % Store ground truth for testing
    params.X0 = X0_ref;
    params.A0 = A1_used;

    for k = 1:params.num_kernels
        params.xinit{k}.X = X_ref(:,:,k);
        params.xinit{k}.b = extras.phase1.biter(k); 
    end


    % Run SBD on reference slice
    fprintf('Finding optimal activation for reference slice %d...\n', params.ref_slice);
    [Aout_slice, Xout_slice, bout_slice, slice_extras] = SBD_test_multi_demixing(...
        Y_used, kernel_sizes_used, params, dispfun, A1_used, initial_iteration, maxIT);

    A_all_out(slice,:) = Aout_slice;
    X_all_out(:,:,:,slice) = Xout_slice;
    b_all_out(:,slice) = bout_slice;
    extras_all_out{slice} = slice_extras;

    fprintf('Completed slice %d/%d\n', slice, params.num_slices);
end 
%%
visualizeResults(Y_used, A1_used, Aout_slice, X0_ref, Xout_slice, bout_slice, extras);


%% Block 5: Run SBD in Parallel for All Slices
% This block runs SBD optimization for all slices in parallel

fprintf('\nRunning SBD optimization for all slices...\n');

num_slices = params.num_slices;
num_kernels = params.num_kernels;
% Pre-process all params_local structures outside parfor
params_all = cell(num_slices, 1);
for n = 1:num_slices
    params_all{n} = struct();
    params_all{n}.lambda1 = params.lambda1;
    params_all{n}.kplus = params.kplus;
    params_all{n}.lambda2 = params.lambda2;
    params_all{n}.xpos = true;  % Non-negative X constraint
    params_all{n}.getbias = true;  % Extract constant bias
    params_all{n}.phase2 = false;
    params_all{n}.Xsolve = 'FISTA';
    params_all{n}.X0 = X0;
    params_all{n}.A0 = A1_all(n,:);
    params_all{n}.xinit = cell(1,num_kernels);
    for k = 1:num_kernels
        sbd_params.xinit{k}.X = X_ref(:,:,k);
        sbd_params.xinit{k}.b = extras.phase1.biter(k); 
    end
end

% Get kernel sizes for all kernels (read-only, safe to share)
kernel_sizes = zeros(params.num_kernels, 2);
for k = 1:params.num_kernels
    kernel_sizes(k,:) = size(A1_all{1,k});
end

% Set up display function (read-only, safe to share)
dispfun = @(Y,A,X,k,kplus,idx) 0;  % No display during parallel processing

% Initialize arrays for results (these will be written to in parallel)
A_all_out = cell(params.num_slices, params.num_kernels);
X_all_out = zeros(size(Y,1), size(Y,2), params.num_kernels, params.num_slices);
b_all_out = zeros(params.num_kernels, params.num_slices);
extras_all_out = cell(params.num_slices, 1);

Y_slices = cell(params.num_slices, 1);
for s = 1:params.num_slices
    Y_slices{s} = Y(:,:,s);
end 

% Run SBD in parallel for all slices
parfor s = num_slices
    fprintf('Processing slice %d/%d...\n', s, num_slices);
    
    % Create local copy of initial kernels for this slice
    A_init_slice = A1_all(s,:);
    % Run SBD for this slice
    [A_final(s,:), X_final(:,:,:,s), b_final(:,s), extras{s}] = ...
        SBD_test_multi_demixing(Y_slices{s}, kernel_sizes, params_all{s}, dispfun, ...
        A_init_slice, initial_iteration, maxIT);
    
    fprintf('Completed slice %d/%d\n', s, num_slices);
end

% Store final results
final_results = struct();
final_results.A = A_final;
final_results.X = X_final;
final_results.b = b_final;
final_results.extras = extras;

fprintf('SBD optimization complete for all slices.\n');

%% Visualize results for reference slice
figure('Name', 'Final Results (Reference Slice)');
for k = 1:params.num_kernels
    % Show final kernel
    subplot(2, params.num_kernels, k);
    imagesc(A_final{params.ref_slice,k});
    colormap(gray);
    colorbar;
    title(sprintf('Final Kernel %d', k));
    axis square;
    
    % Show final activation map
    subplot(2, params.num_kernels, k + params.num_kernels);
    imagesc(X_final(:,:,k,params.ref_slice));
    colormap(gray);
    colorbar;
    title(sprintf('Final Activation %d', k));
    axis square;
end

%% Block 6: Analyze Results
% This block analyzes the results of the SBD optimization

fprintf('\nAnalyzing results...\n');

% Calculate average activation values for each kernel
average_activations = zeros(params.num_kernels, 1);
for k = 1:params.num_kernels
    average_activations(k) = mean(X_final(:,:,k,params.ref_slice), 'all');
end

% Calculate average kernel sizes for each kernel
average_kernel_sizes = zeros(params.num_kernels, 2);
for k = 1:params.num_kernels
    average_kernel_sizes(k,:) = mean(kernel_sizes(k,:), 2);
end

% Calculate average isolation scores for each kernel
average_isolation_scores = zeros(params.num_kernels, 1);
for k = 1:params.num_kernels
    average_isolation_scores(k) = mean(isolation_scores{k});
end

% Display results
fprintf('Average activation values for each kernel:\n');
disp(average_activations);

fprintf('Average kernel sizes for each kernel:\n');
disp(average_kernel_sizes);

fprintf('Average isolation scores for each kernel:\n');
disp(average_isolation_scores);

fprintf('Analysis complete.\n');

%% Helper Functions
% These functions are used by the optimization blocks

function [cost, store] = costfun_proximity(a, store, Y, X_opt, lambda, mu, lambda_prox)
    % Cost function with proximity regularization
    
    % Reshape a to kernel
    A = reshape(a, size(X_opt));
    
    % Compute X using FISTA
    if ~isfield(store, 'X')
        [store.X, ~] = Xsolve_FISTA_tunable(Y, A, lambda, mu);
    end
    
    % Data fidelity term
    data_fit = norm(convfft2(A, store.X) - Y, 'fro')^2/2;
    
    % Regularization term (using Huber function)
    reg_term = huber_cost(store.X, mu, lambda);
    
    % Proximity term
    prox_term = norm(store.X - X_opt, 'fro')^2;
    
    % Combined cost
    cost = data_fit + reg_term + lambda_prox * prox_term;
end

function [egrad, store] = egradfun_proximity(a, store, Y, X_opt, lambda, mu, lambda_prox)
    % Euclidean gradient with proximity regularization
    
    % Reshape a to kernel
    A = reshape(a, size(X_opt));
    
    % Ensure X is computed
    if ~isfield(store, 'X')
        [store.X, ~] = Xsolve_FISTA_tunable(Y, A, lambda, mu);
    end
    
    % Data fidelity gradient
    data_grad = convfft2(store.X, convfft2(A, store.X) - Y);
    
    % Proximity gradient (through X)
    [~, prox_grad_X] = huber_grad(store.X, mu, lambda);
    prox_grad = 2 * lambda_prox * (store.X - X_opt);
    
    % Combined gradient
    egrad = data_grad(:) + prox_grad(:);
end

function [ehess, store] = ehessfun_proximity(a, u, store, Y, X_opt, lambda, mu, lambda_prox)
    % Hessian approximation with proximity regularization
    
    % Reshape vectors to matrices
    A = reshape(a, size(X_opt));
    U = reshape(u, size(X_opt));
    
    % Ensure X is computed
    if ~isfield(store, 'X')
        [store.X, ~] = Xsolve_FISTA_tunable(Y, A, lambda, mu);
    end
    
    % Data fidelity Hessian
    hess_data = convfft2(store.X, convfft2(U, store.X));
    
    % Proximity Hessian
    hess_prox = 2 * lambda_prox * U;
    
    % Combined Hessian
    ehess = hess_data(:) + hess_prox(:);
end

function cost = huber_cost(X, mu, lambda)
    % Huber cost function
    abs_X = abs(X);
    cost = sum(sum(lambda * (abs_X.^2/2) .* (abs_X <= mu) + ...
                  lambda * (mu * abs_X - mu^2/2) .* (abs_X > mu)));
end

function [grad, hess] = huber_grad(X, mu, lambda)
    % Huber gradient and Hessian
    abs_X = abs(X);
    grad = lambda * X .* (abs_X <= mu) + ...
           lambda * mu * sign(X) .* (abs_X > mu);
    
    if nargout > 1
        hess = lambda * (abs_X <= mu);
    end
end 