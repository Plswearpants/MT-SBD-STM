%% ========================================================================
%  TRUNK SCRIPT: Synthetic Data MT-SBD-STM Analysis
%  ========================================================================
%  Multi-kernel Tensor Shifted Blind Deconvolution for Scanning Tunneling 
%  Microscopy synthetic data generation, processing, and visualization.
%
%  This script follows the UBC LAIR standardized block structure for
%  reproducibility and comprehensive logging.
%
%  Author: Dong Chen
%  Created: 2025-10-27
%  Project: MT-SBD-STM
%
%  WORKFLOW OVERVIEW:
%  ==================
%  Block 1 (GD01A): Generate synthetic test data
%  Block 2 (IN01A): Initialize kernels for reference slice
%  Block 3 (DS01A): Decompose reference slice (find optimal activation)
%  Block 4 (IS01A): Find most isolated points for kernel proliferation
%  Block 5 (IP01A): Initialize kernels for all slices (proliferation)
%  Block 6 (DA01A): Decompose all slices simultaneously
%  Block 7 (VR01A): Visualize results
%  Block 8 (WS01A): Save workspace and results
%
%  BLOCK ID CONVENTION: ABXXZ
%  ==========================
%  A = Category:    G=Generate, I=Initialize, D=Decompose, V=Visualize, 
%                   W=Write, S=Select, P=Process
%  B = Subcategory: D=Data, N=kNernel, S=Slice, P=Proliferation, 
%                   R=Results, W=Workspace, A=Analysis
%  XX = Number:     01, 02, 03, ...
%  Z = Variant:     A, B, C, ...
%
%  PARAMETER GLOSSARY:
%  ===================
%  See docs/parameter_glossary.md for complete parameter definitions
%
% =========================================================================

%% SECTION 0: WORKSPACE INITIALIZATION
% =========================================================================
% Clear workspace and initialize paths
clc; clear; close all;

% Add paths (run init_sbd from parent directory)
if exist('../init_sbd.m', 'file')
    run('../init_sbd.m');
else
    error('init_sbd.m not found. Please run from scripts/ directory.');
end

% Initialize logging system
LOGpath = pwd;  % Current directory for logs
LOGfile = sprintf('synthetic_run_%s', datestr(now, 'yyyymmdd_HHMMSS'));
LOGcomment = "";

% Initialize log file (clear if exists)
[LOGcomment] = logUsedBlocks(LOGpath, LOGfile, "INIT ", "=== Synthetic Data Workflow Started ===", 1);

fprintf('\n========================================\n');
fprintf('MT-SBD-STM Synthetic Data Analysis\n');
fprintf('========================================\n');
fprintf('Log file: %s_LOGfile.txt\n\n', LOGfile);


%% =========================================================================
%% GD01A: Generate-Data-01-A; Generate Synthetic Test Data
%  =========================================================================
%  Edited by Dong Chen, 2025-10-27
%
%  Generates synthetic 3D STM observation data with multiple kernels across
%  energy slices. Uses LDoS simulation data and creates realistic defect
%  patterns with controlled SNR.
%
%  Dependencies: generateSyntheticData.m (wrapper)

% -------------------------------------------------------------------------
% PRESETS: User-configurable parameters
% -------------------------------------------------------------------------
SNR = 5;                        % Signal-to-noise ratio
N_obs = 50;                     % Observation lattice size (pixels)
observation_resolution = 3;     % Resolution: pixels per lattice site
defect_density = 0.01;          % Surface defect density (0-1)
num_slices = 2;                 % Number of energy slices
LDoS_path = 'LDoS_single_defects_self=0.6_save.mat';
vis_generation = false;         % Show intermediate generation steps
normalization_type = 'dynamic'; % 'dynamic' or 'static'
ref_slice = [];                 % Reference slice ([] for interactive selection)

%%%%%%%%%%%%%%%%%% DO NOT EDIT BELOW %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% PRESET LOG: 
LOGcomment = sprintf("SNR=%g, N_obs=%d, resolution=%d, density=%g, slices=%d", ...
    SNR, N_obs, observation_resolution, defect_density, num_slices);
LOGcomment = logUsedBlocks(LOGpath, LOGfile, "GD01A", LOGcomment, 0);

% Generate synthetic data (compact 1-line function call)
fprintf('Generating synthetic test data...\n');
[data, params] = generateSyntheticData(...
    'SNR', SNR, 'N_obs', N_obs, 'observation_resolution', observation_resolution, ...
    'defect_density', defect_density, 'num_slices', num_slices, 'LDoS_path', LDoS_path, ...
    'vis_generation', vis_generation, 'normalization_type', normalization_type, ...
    'ref_slice', ref_slice);

fprintf('Data generation complete.\n\n');

% LOG: Generation details and reference slice
LOGcomment = sprintf("Generated Y: %dx%dx%d, Kernels: %d, Ref slice: %d", ...
    size(data.Y,1), size(data.Y,2), size(data.Y,3), params.num_kernels, params.ref_slice);
LOGcomment = logUsedBlocks(LOGpath, LOGfile, "  ^  ", LOGcomment, 0);

% Clear preset variables
clearvars SNR N_obs observation_resolution defect_density num_slices LDoS_path
clearvars vis_generation normalization_type ref_slice


%% =========================================================================
%% IN01A: Initialize-kNernel-01-A; Initialize Kernels for Reference Slice
%  =========================================================================
%  Edited by Dong Chen, 2025-10-27
%
%  Interactive kernel initialization from the reference slice observation.
%  User selects regions containing isolated defect patterns to use as
%  initial kernel guesses.
%
%  Dependencies: initialize_kernels.m, apply_window.m

% -------------------------------------------------------------------------
% PRESETS: User-configurable parameters
% -------------------------------------------------------------------------
kernel_selection_type = 'selected';     % 'selected' or 'random'
window_type = {'gaussian', 2.5};        % Window function for kernels
                                        % Options: 'hann', 'hamming', 'blackman'
                                        %          {'gaussian', alpha}
                                        %          {'kaiser', beta}
                                        %          {} for no window

%%%%%%%%%%%%%%%%%% DO NOT EDIT BELOW %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% LOG: Block start
LOGcomment = sprintf("Selection: %s, Window: %s", kernel_selection_type, mat2str(window_type));
LOGcomment = logUsedBlocks(LOGpath, LOGfile, "IN01A", LOGcomment, 0);

% Get kernel sizes for reference slice (temporary for this block)
kernel_sizes_ref_temp = reshape(params.kernel_sizes(params.ref_slice,:,:), [params.num_kernels, 2]);

% Display ground truth kernels
figure('Name', 'IN01A: Ground Truth Kernels');
for k = 1:params.num_kernels
    subplot(1, params.num_kernels, k);
    imagesc(A0{params.ref_slice, k});
    axis square;
    title(sprintf('True Kernel %d', k));
    colorbar;
end
sgtitle('Ground Truth Kernels (Reference Slice)');

% Initialize kernels interactively
fprintf('Initializing kernels for reference slice...\n');
fprintf('Please select %d kernel regions from the observation.\n', params.num_kernels);
A1 = initialize_kernels(Y_ref, params.num_kernels, kernel_sizes_ref_temp, kernel_selection_type, window_type);

% Display initialized kernels
figure('Name', 'IN01A: Initialized Kernels');
for n = 1:params.num_kernels
    subplot(1, params.num_kernels, n);
    imagesc(A1{n});
    title(sprintf('Initial Kernel %d', n));
    colorbar;
    axis square;
end
sgtitle('Initialized Kernels (Reference Slice)');

fprintf('Kernel initialization complete.\n\n');

% LOG: Initialization complete
LOGcomment = sprintf("Initialized %d kernels, sizes: %s", params.num_kernels, mat2str(kernel_sizes_ref_temp));
LOGcomment = logUsedBlocks(LOGpath, LOGfile, "  ^  ", LOGcomment, 0);

% Clear preset variables
clearvars kernel_selection_type window_type kernel_sizes_ref_temp n k


%% =========================================================================
%% DS01A: Decompose-Slice-01-A; Decompose Reference Slice
%  =========================================================================
%  Edited by Dong Chen, 2025-10-27
%
%  Runs the MT-SBD algorithm on the reference slice to find optimal
%  activation maps and refine kernel estimates. This provides the starting
%  point for multi-slice analysis.
%
%  Dependencies: MTSBD_synthetic.m, Asolve_Manopt_tunable.m, 
%                Xsolve_FISTA_tunable.m, computeQualityMetrics.m, showims.m

% -------------------------------------------------------------------------
% PRESETS: User-configurable parameters
% -------------------------------------------------------------------------
% Algorithm selection
use_Xregulated = false;         % true: MTSBD_Xregulated, false: MTSBD_synthetic

% Phase I settings
initial_iteration = 1;          % Manopt/FISTA inner iterations (start)
maxIT = 15;                     % Number of outer alternating iterations
lambda1 = [3e-2, 3e-2, 3e-2, 3e-2];  % L1 regularization (Phase I)

% Phase II settings (refinement)
phase2_enable = false;          % Enable Phase II refinement
lambda2 = [1e-2, 1e-2, 1e-2, 1e-2];  % L1 regularization (Phase II final)
nrefine = 5;                    % Number of refinement steps
kplus_factor = 0.5;             % Sphere lifting padding factor

% Algorithm parameters
signflip_threshold = 0.2;       % Sign flip detection threshold
xpos = true;                    % Enforce positive activations
getbias = true;                 % Extract constant bias term
Xsolve_method = 'FISTA';        % 'FISTA' or 'pdNCG'
gamma_crosscorr = 5e-2;         % Cross-correlation regularization

% Initialization options
use_xinit = [];                 % Initial X guess ([] for none)

% Display options
show_progress = true;           % Show progress during optimization

%%%%%%%%%%%%%%%%%% DO NOT EDIT BELOW %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% LOG: Block start
LOGcomment = sprintf("Method: %s, maxIT=%d, lambda1=%s, phase2=%d", ...
    iif(use_Xregulated, 'Xregulated', 'standard'), maxIT, mat2str(lambda1), phase2_enable);
LOGcomment = logUsedBlocks(LOGpath, LOGfile, "DS01A", LOGcomment, 0);

% Set up display functions for monitoring
if show_progress
    figure('Name', 'DS01A: Decomposition Progress');
    dispfun = cell(1, num_kernels);
    for n = 1:num_kernels
        dispfun{n} = @(Y, A, X, kernel_sizes, kplus) showims(Y_ref, A1{n}, X0_ref(:,:,n), A, X, kernel_sizes, kplus, 1);
    end
else
    dispfun = cell(1, num_kernels);
    for n = 1:num_kernels
        dispfun{n} = @(Y, A, X, kernel_sizes, kplus) 0;
    end
end

% Package parameters
kernel_sizes_ref = reshape(kernel_sizes(params.ref_slice,:,:), [num_kernels, 2]);
sbd_params = struct();
sbd_params.lambda1 = lambda1;
sbd_params.phase2 = phase2_enable;
sbd_params.kplus = ceil(kplus_factor * kernel_sizes_ref);
sbd_params.lambda2 = lambda2;
sbd_params.nrefine = nrefine;
sbd_params.signflip = signflip_threshold;
sbd_params.xpos = xpos;
sbd_params.getbias = getbias;
sbd_params.Xsolve = Xsolve_method;
sbd_params.X0 = X0_ref;
sbd_params.A0 = A0_ref;
sbd_params.xinit = use_xinit;
sbd_params.gamma = gamma_crosscorr;

% Run SBD on reference slice
fprintf('Running MT-SBD decomposition on reference slice %d...\n', params.ref_slice);
fprintf('Method: %s\n', iif(use_Xregulated, 'X-regularized', 'Standard'));
fprintf('Max iterations: %d\n', maxIT);
tic;

if use_Xregulated
    [A_ref, X_ref, bout, extras] = MTSBD_synthetic_Xregulated(...
        Y_ref, kernel_sizes_ref, sbd_params, dispfun, A1, initial_iteration, maxIT);
else
    [A_ref, X_ref, bout, extras] = MTSBD_synthetic(...
        Y_ref, kernel_sizes_ref, sbd_params, dispfun, A1, initial_iteration, maxIT);
end

decomp_time = toc;
fprintf('Decomposition completed in %.2f seconds.\n', decomp_time);

% Store reference results
ref_results = struct();
ref_results.A = A_ref;
ref_results.X = X_ref;
ref_results.b = bout;
ref_results.extras = extras;
ref_results.params = sbd_params;

% Display final quality metrics
fprintf('\nFinal Quality Metrics:\n');
final_metrics = extras.phase1.activation_metrics(end,:);
final_kernel_quality = extras.phase1.kernel_quality_factors(end,:);
for k = 1:num_kernels
    fprintf('  Kernel %d - Activation: %.4f, Quality: %.4f\n', ...
        k, final_metrics(k), final_kernel_quality(k));
end
fprintf('\n');

% LOG: Decomposition results
LOGcomment = sprintf("Completed in %.2fs, Final metrics: Act=%s, Qual=%s", ...
    decomp_time, mat2str(final_metrics, 3), mat2str(final_kernel_quality, 3));
LOGcomment = logUsedBlocks(LOGpath, LOGfile, "  ^  ", LOGcomment, 0);

% Clear preset variables
clearvars use_Xregulated initial_iteration maxIT lambda1 phase2_enable lambda2
clearvars nrefine kplus_factor signflip_threshold xpos getbias Xsolve_method
clearvars gamma_crosscorr use_xinit show_progress dispfun sbd_params
clearvars kernel_sizes_ref decomp_time final_metrics final_kernel_quality k n


%% =========================================================================
%% IS01A: Isolation-Selection-01-A; Find Most Isolated Points
%  =========================================================================
%  Edited by Dong Chen, 2025-10-27
%
%  Analyzes the activation maps from the reference slice to find the most
%  isolated defect for each kernel. These points are used as centers for
%  kernel proliferation across all energy slices.
%
%  Dependencies: alignActivationMaps.m, padKernels.m

% -------------------------------------------------------------------------
% PRESETS: User-configurable parameters
% -------------------------------------------------------------------------
isolation_threshold_factor = 10;    % Threshold = max/factor for defect detection
target_kernel_size_type = 'kernel_sizes_all';  % 'ref_kernel_sizes', 
                                               % 'kernel_sizes_cap', 
                                               % 'kernel_sizes_all'

%%%%%%%%%%%%%%%%%% DO NOT EDIT BELOW %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% LOG: Block start
LOGcomment = sprintf("Threshold factor=%g, Size type=%s", ...
    isolation_threshold_factor, target_kernel_size_type);
LOGcomment = logUsedBlocks(LOGpath, LOGfile, "IS01A", LOGcomment, 0);

fprintf('Calculating isolation scores and finding most isolated points...\n');

% Determine target kernel sizes
switch target_kernel_size_type
    case 'ref_kernel_sizes'
        target_kernel_sizes = squeeze(kernel_sizes(params.ref_slice,:,:));
    case 'kernel_sizes_cap'
        target_kernel_sizes = squeeze(max(kernel_sizes,[],1));
    case 'kernel_sizes_all'
        target_kernel_sizes = kernel_sizes;
end

% Initialize storage
most_isolated_points = cell(1, num_kernels);
isolation_scores = cell(1, num_kernels);
defect_positions = cell(1, num_kernels);
num_defects = zeros(1, num_kernels);

% Analyze activation value distributions
figure('Name', 'IS01A: Activation Value Distributions');
for k = 1:num_kernels
    % Histogram of activation values
    subplot(2, num_kernels, k);
    activation_values = X_ref(:,:,k);
    histogram(activation_values(activation_values > 0), 50);
    set(gca, 'YScale', 'log');
    title(sprintf('Kernel %d Distribution', k));
    xlabel('Activation Value');
    ylabel('Frequency (log)');
    
    % Threshold line
    threshold = max(X_ref(:,:,k), [], 'all') / isolation_threshold_factor;
    hold on;
    xline(threshold, 'r--', 'Threshold');
    hold off;
    
    % Cumulative distribution
    subplot(2, num_kernels, k + num_kernels);
    [counts, edges] = histcounts(activation_values(activation_values > 0), 50, 'Normalization', 'cdf');
    stairs(edges(1:end-1), counts);
    title(sprintf('Kernel %d CDF', k));
    xlabel('Activation Value');
    ylabel('Cumulative Frequency');
    hold on;
    xline(threshold, 'r--', 'Threshold');
    hold off;
    
    % Get defect positions above threshold
    [rows, cols] = find(X_ref(:,:,k) > threshold);
    defect_positions{k} = [rows, cols];
    num_defects(k) = size(defect_positions{k}, 1);
    fprintf('  Kernel %d: %d defects above threshold %.4f\n', k, num_defects(k), threshold);
end

% Calculate isolation scores
for k = 1:num_kernels
    if num_defects(k) == 0
        warning('No defects found for kernel %d', k);
        continue;
    end
    
    % Sum activation maps of all other kernels
    X_others = zeros(size(X_ref(:,:,1)));
    for l = 1:num_kernels
        if l ~= k
            X_others = X_others + X_ref(:,:,l);
        end
    end
    
    % Get positions in other kernels
    [other_rows, other_cols] = find(X_others > max(X_others,[],'all') / isolation_threshold_factor);
    other_positions = [other_rows, other_cols];
    
    % Boundary check based on target kernel size
    if strcmp(target_kernel_size_type, 'kernel_sizes_all')
        half_kernel_size = floor(squeeze(target_kernel_sizes(params.ref_slice,k,:))' / 2);
    else
        half_kernel_size = floor(target_kernel_sizes(k,:) / 2);
    end
    
    % Filter out points too close to boundaries
    valid_points = true(num_defects(k), 1);
    for i = 1:num_defects(k)
        y = defect_positions{k}(i,1);
        x = defect_positions{k}(i,2);
        
        if y <= half_kernel_size(1) || y >= size(X_ref,1) - half_kernel_size(1) || ...
           x <= half_kernel_size(2) || x >= size(X_ref,2) - half_kernel_size(2)
            valid_points(i) = false;
        end
    end
    
    % Process valid points
    valid_defects = defect_positions{k}(valid_points,:);
    if isempty(valid_defects)
        error('No valid isolated points for kernel %d - all too close to boundaries', k);
    end
    
    % Calculate isolation scores (distance to nearest other-kernel defect)
    S_k = zeros(size(valid_defects, 1), 1);
    for i = 1:size(valid_defects, 1)
        diffs = other_positions - valid_defects(i,:);
        distances = sum(diffs.^2, 2);
        S_k(i) = min(distances);
    end
    
    % Find most isolated point
    [max_score, max_idx] = max(S_k);
    valid_indices = find(valid_points);
    max_idx = valid_indices(max_idx);
    
    most_isolated_points{k} = defect_positions{k}(max_idx,:);
    isolation_scores{k} = S_k;
    
    fprintf('  Kernel %d: Most isolated at (%d,%d), score=%.2f\n', ...
        k, most_isolated_points{k}(1), most_isolated_points{k}(2), max_score);
end

% Visualize isolation analysis
figure('Name', 'IS01A: Isolation Analysis');
for k = 1:num_kernels
    % Activation map with all defects and most isolated point
    subplot(2, num_kernels, k);
    imagesc(X_ref(:,:,k));
    hold on;
    scatter(defect_positions{k}(:,2), defect_positions{k}(:,1), 50, 'w', 'o');
    if ~isempty(most_isolated_points{k})
        scatter(most_isolated_points{k}(2), most_isolated_points{k}(1), 100, 'r', '*', 'LineWidth', 2);
    end
    title(sprintf('Kernel %d', k));
    colorbar;
    axis square;
    hold off;
    
    % Other kernels with most isolated point marked
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
        scatter(most_isolated_points{k}(2), most_isolated_points{k}(1), 100, 'r', '*', 'LineWidth', 2);
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
isolation_analysis.target_kernel_sizes = target_kernel_sizes;

% Prepare ground truth kernels with appropriate sizes
if strcmp(target_kernel_size_type, 'kernel_sizes_all')
    A0_used = A0;
else
    A0_used = padKernels(A0_noiseless, params.SNR, target_kernel_sizes);
end

% Align most isolated points with ground truth
kernel_sizes_ref_temp = reshape(kernel_sizes(params.ref_slice,:,:), [num_kernels, 2]);
[~, offset, ~] = alignActivationMaps(X0_ref, X_ref, kernel_sizes_ref_temp);
used_most_isolated_points = cell(1, num_kernels);
for k = 1:num_kernels
    used_most_isolated_points{k} = most_isolated_points{k} + offset(k,:);
end

% Display most isolated points on observation
figure('Name', 'IS01A: Isolated Points on Observation');
imagesc(Y_ref);
colormap(gray);
hold on;
for k = 1:num_kernels
    scatter(most_isolated_points{k}(2), most_isolated_points{k}(1), 100, 'r', '*', 'LineWidth', 2);
    text(most_isolated_points{k}(2)+5, most_isolated_points{k}(1), sprintf('K%d', k), ...
        'Color', 'red', 'FontWeight', 'bold', 'FontSize', 12);
end
hold off;
title('Most Isolated Points for Each Kernel');
axis square;
colorbar;

fprintf('Isolation analysis complete.\n\n');

% LOG: Isolation results
LOGcomment = sprintf("Found %d isolated points, offsets: %s", ...
    num_kernels, mat2str(offset));
LOGcomment = logUsedBlocks(LOGpath, LOGfile, "  ^  ", LOGcomment, 0);

% Clear preset variables
clearvars isolation_threshold_factor target_kernel_size_type threshold
clearvars activation_values counts edges rows cols other_rows other_cols
clearvars half_kernel_size valid_points valid_defects S_k max_score max_idx
clearvars valid_indices X_others diffs distances i k l
clearvars kernel_sizes_ref_temp offset


%% =========================================================================
%% IP01A: Initialize-Proliferation-01-A; Initialize Kernels for All Slices
%  =========================================================================
%  Edited by Dong Chen, 2025-10-27
%
%  Uses the most isolated points identified in the reference slice as
%  centers to initialize kernels across all energy slices. This "proliferation"
%  approach ensures consistent spatial positioning across the energy dimension.
%
%  Dependencies: initialize_kernels_proliferation.m

% -------------------------------------------------------------------------
% PRESETS: User-configurable parameters
% -------------------------------------------------------------------------
window_type_proliferation = {'gaussian', 2.5};  % Window function
interactive_size_adjust = false;                % Allow interactive resizing
use_matrix_format = true;                       % Output as matrix vs cell array

%%%%%%%%%%%%%%%%%% DO NOT EDIT BELOW %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% LOG: Block start
LOGcomment = sprintf("Window: %s, Interactive: %d, Matrix: %d", ...
    mat2str(window_type_proliferation), interactive_size_adjust, use_matrix_format);
LOGcomment = logUsedBlocks(LOGpath, LOGfile, "IP01A", LOGcomment, 0);

fprintf('Initializing kernels for all %d slices using proliferation...\n', num_slices);

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

for s = 1:num_slices
    fprintf('  Slice %d/%d...', s, num_slices);
    
    if strcmp(isolation_analysis.target_kernel_sizes, 'kernel_sizes_all') || ndims(target_kernel_sizes) == 3
        % Use slice-specific kernel sizes
        target_sizes_slice = squeeze(kernel_sizes(s,:,:));
    else
        % Use uniform kernel sizes
        target_sizes_slice = target_kernel_sizes;
    end
    
    A1_all(s,:) = initialize_kernels_proliferation(Y(:,:,s), num_kernels, ...
        kernel_centers, window_type_proliferation, target_sizes_slice, ...
        'interactive', interactive_size_adjust);
    
    fprintf(' done.\n');
end

% Visualize initialized kernels for reference slice
figure('Name', 'IP01A: Initialized Kernels (Reference Slice)');
for k = 1:num_kernels
    subplot(1, num_kernels, k);
    imagesc(A1_all{params.ref_slice, k});
    colormap(gray);
    colorbar;
    title(sprintf('Initialized Kernel %d', k));
    axis square;
end
sgtitle(sprintf('Proliferated Kernels - Slice %d', params.ref_slice));

% Store initialization results
init_results = struct();
init_results.A = A1_all;
init_results.kernel_centers = kernel_centers;
init_results.kernel_sizes = target_kernel_sizes;

fprintf('Kernel initialization complete for all slices.\n\n');

% Convert to matrix format if requested
if use_matrix_format
    A1_all_matrix = cell(num_kernels, 1);
    for k = 1:num_kernels
        A1_all_matrix{k} = zeros(size(A1_all{1,k},1), size(A1_all{1,k},2), num_slices);
        for s = 1:num_slices
            A1_all_matrix{k}(:,:,s) = A1_all{s,k};
        end
    end
    fprintf('Converted to matrix format: %d kernels x [%dx%dx%d]\n\n', ...
        num_kernels, size(A1_all_matrix{1},1), size(A1_all_matrix{1},2), num_slices);
end

% LOG: Initialization complete
LOGcomment = sprintf("Initialized %d kernels for %d slices at centers: %s", ...
    num_kernels, num_slices, mat2str(kernel_centers));
LOGcomment = logUsedBlocks(LOGpath, LOGfile, "  ^  ", LOGcomment, 0);

% Clear preset variables
clearvars window_type_proliferation interactive_size_adjust use_matrix_format
clearvars kernel_centers target_sizes_slice s k


%% =========================================================================
%% DA01A: Decompose-All-01-A; Decompose All Slices Simultaneously
%  =========================================================================
%  Edited by Dong Chen, 2025-10-27
%
%  Runs MT-SBD on all energy slices simultaneously, using the reference
%  slice activation as initialization. This leverages 3D convolution for
%  efficient multi-slice processing.
%
%  Dependencies: MTSBD_synthetic_all_slice.m, convfft3.m

% -------------------------------------------------------------------------
% PRESETS: User-configurable parameters
% -------------------------------------------------------------------------
% Algorithm selection
use_Xregulated_allslice = false;    % Use X-regularized version

% Use reference slice results as initialization
use_reference_init = true;          % Initialize X with reference results

% Display options
show_allslice_progress = true;      % Show progress during optimization

% Algorithm parameters inherited from DS01A
maxIT_allslice = 15;                % Max iterations for all-slice decomposition

%%%%%%%%%%%%%%%%%% DO NOT EDIT BELOW %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% LOG: Block start
LOGcomment = sprintf("Method: %s, Init from ref: %d, maxIT=%d", ...
    iif(use_Xregulated_allslice, 'Xregulated', 'standard'), use_reference_init, maxIT_allslice);
LOGcomment = logUsedBlocks(LOGpath, LOGfile, "DA01A", LOGcomment, 0);

fprintf('Running MT-SBD decomposition on all %d slices...\n', num_slices);

% Prepare kernel sizes and observation
kernel_sizes_used = squeeze(max(kernel_sizes, [], 1));
Y_used = Y;
A1_used = A1_all_matrix;

% Set up display functions
if show_allslice_progress
    figure('Name', 'DA01A: All-Slice Decomposition Progress');
    dispfun_allslice = cell(1, num_kernels);
    for n = 1:num_kernels
        dispfun_allslice{n} = @(Y, A, X, kernel_sizes, kplus) ...
            showims(Y_used, A1_used{n}, X0_ref(:,:,n), A, X, kernel_sizes, kplus, 1);
    end
else
    dispfun_allslice = cell(1, num_kernels);
    for n = 1:num_kernels
        dispfun_allslice{n} = @(Y, A, X, kernel_sizes, kplus) 0;
    end
end

% Update parameters for all-slice processing
params_allslice = ref_results.params;
params_allslice.X0 = X0;
params_allslice.A0 = A0_used;

% Initialize X from reference slice if requested
if use_reference_init
    params_allslice.xinit = cell(1, num_kernels);
    for k = 1:num_kernels
        params_allslice.xinit{k}.X = X_ref(:,:,k);
        b_temp = ref_results.extras.phase1.biter(k);
        params_allslice.xinit{k}.b = repmat(b_temp, [num_slices, 1]);
    end
else
    params_allslice.xinit = [];
end

% Run all-slice decomposition
tic;
if use_Xregulated_allslice
    [Aout_slice, Xout_slice, bout_slice, slice_extras] = MTSBD_synthetic_Xregulated_all_slices(...
        Y_used, kernel_sizes_used, params_allslice, dispfun_allslice, A1_used, initial_iteration, maxIT_allslice);
else
    [Aout_slice, Xout_slice, bout_slice, slice_extras] = MTSBD_synthetic_all_slice(...
        Y_used, kernel_sizes_used, params_allslice, dispfun_allslice, A1_used, initial_iteration, maxIT_allslice);
end
allslice_time = toc;

fprintf('All-slice decomposition completed in %.2f seconds.\n', allslice_time);

% Store all-slice results
allslice_results = struct();
allslice_results.A = Aout_slice;
allslice_results.X = Xout_slice;
allslice_results.b = bout_slice;
allslice_results.extras = slice_extras;
allslice_results.params = params_allslice;

fprintf('Multi-slice decomposition complete.\n\n');

% LOG: All-slice results
LOGcomment = sprintf("Completed in %.2fs for %d slices", allslice_time, num_slices);
LOGcomment = logUsedBlocks(LOGpath, LOGfile, "  ^  ", LOGcomment, 0);

% Clear preset variables
clearvars use_Xregulated_allslice use_reference_init show_allslice_progress
clearvars maxIT_allslice kernel_sizes_used Y_used A1_used dispfun_allslice
clearvars params_allslice b_temp allslice_time n k


%% =========================================================================
%% VR01A: Visualize-Results-01-A; Visualize All Results
%  =========================================================================
%  Edited by Dong Chen, 2025-10-27
%
%  Comprehensive visualization of decomposition results for all slices,
%  including reconstructions, activations, kernels, and quality metrics.
%
%  Dependencies: visualizeResults.m

% -------------------------------------------------------------------------
% PRESETS: User-configurable parameters
% -------------------------------------------------------------------------
visualize_all_slices = true;        % Visualize every slice (false for ref only)
save_figures = false;               % Save figure files

%%%%%%%%%%%%%%%%%% DO NOT EDIT BELOW %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% LOG: Block start
LOGcomment = sprintf("Visualize all: %d, Save: %d", visualize_all_slices, save_figures);
LOGcomment = logUsedBlocks(LOGpath, LOGfile, "VR01A", LOGcomment, 0);

fprintf('Generating result visualizations...\n');

% Convert Aout_slice to cell format for visualization
Aout_cell = cell(num_slices, num_kernels);
for s = 1:num_slices
    for k = 1:num_kernels
        Aout_cell{s,k} = Aout_slice{k}(:,:,s);
    end
end

% Determine which slices to visualize
if visualize_all_slices
    slices_to_viz = 1:num_slices;
else
    slices_to_viz = params.ref_slice;
end

% Visualize results for each selected slice
for s = slices_to_viz
    fprintf('  Visualizing slice %d/%d...\n', s, num_slices);
    
    % Prepare slice-specific data
    Y_slice = Y(:,:,s);
    A0_slice = cell(1, num_kernels);
    for k = 1:num_kernels
        if iscell(A0_used)
            A0_slice{k} = A0_used{k}(:,:,s);
        else
            A0_slice{k} = A0_used{s,k};
        end
    end
    Aout_slice_cell = cell(1, num_kernels);
    for k = 1:num_kernels
        Aout_slice_cell{k} = Aout_cell{s,k};
    end
    
    % Call visualizeResults for this slice
    visualizeResults(Y_slice, A0_slice, Aout_slice_cell, X0, Xout_slice, ...
        bout_slice(s,:), slice_extras, [s, 1]);
    
    % Save figures if requested
    if save_figures
        fig_name = sprintf('%s_slice%02d_results', LOGfile, s);
        savefig(gcf, fullfile(LOGpath, [fig_name '.fig']));
        LOGcomment = sprintf("Saved figure: %s.fig", fig_name);
        LOGcomment = logUsedBlocks(LOGpath, LOGfile, "  ^  ", LOGcomment, 0);
    end
end

fprintf('Visualization complete.\n\n');

% LOG: Visualization complete
LOGcomment = sprintf("Visualized %d slices", length(slices_to_viz));
LOGcomment = logUsedBlocks(LOGpath, LOGfile, "  ^  ", LOGcomment, 0);

% Clear preset variables
clearvars visualize_all_slices save_figures slices_to_viz fig_name
clearvars Y_slice A0_slice Aout_slice_cell s k


%% =========================================================================
%% WS01A: Write-Save-01-A; Save Workspace and Results
%  =========================================================================
%  Edited by Dong Chen, 2025-10-27
%
%  Saves the complete workspace including all results, parameters, and
%  intermediate data for later analysis. Also saves the log file.

% -------------------------------------------------------------------------
% PRESETS: User-configurable parameters
% -------------------------------------------------------------------------
save_workspace = true;              % Save complete workspace
save_results_only = false;          % Save only key results (smaller file)
results_path = pwd;                 % Path for saving results

%%%%%%%%%%%%%%%%%% DO NOT EDIT BELOW %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% LOG: Block start
LOGcomment = sprintf("Save workspace: %d, Results only: %d", save_workspace, save_results_only);
LOGcomment = logUsedBlocks(LOGpath, LOGfile, "WS01A", LOGcomment, 0);

if save_workspace || save_results_only
    fprintf('Saving results...\n');
    
    % Create results filename
    results_filename = sprintf('%s_results.mat', LOGfile);
    results_fullpath = fullfile(results_path, results_filename);
    
    if save_results_only
        % Save only essential results
        save(results_fullpath, 'Y', 'A0', 'X0', 'params', ...
            'ref_results', 'allslice_results', 'isolation_analysis', ...
            'num_kernels', 'num_slices', 'LOGfile', '-v7.3');
        fprintf('Results saved to: %s\n', results_fullpath);
    else
        % Save complete workspace
        save(results_fullpath, '-v7.3');
        fprintf('Complete workspace saved to: %s\n', results_fullpath);
    end
    
    % Copy log file to results location
    log_source = fullfile(LOGpath, [LOGfile '_LOGfile.txt']);
    log_dest = fullfile(results_path, [LOGfile '_LOGfile.txt']);
    if ~strcmp(log_source, log_dest)
        copyfile(log_source, log_dest);
    end
    
    fprintf('Log file saved to: %s\n', log_dest);
    fprintf('\n');
    
    % LOG: Save complete
    LOGcomment = sprintf("Saved to: %s", results_filename);
    LOGcomment = logUsedBlocks(LOGpath, LOGfile, "  ^  ", LOGcomment, 0);
end

% Clear preset variables
clearvars save_workspace save_results_only results_path
clearvars results_filename results_fullpath log_source log_dest


%% =========================================================================
%% FINAL: Analysis Complete
%  =========================================================================

% LOG: Final entry
LOGcomment = "=== Synthetic Data Analysis Complete ===";
LOGcomment = logUsedBlocks(LOGpath, LOGfile, "DONE ", LOGcomment, 0);

fprintf('========================================\n');
fprintf('MT-SBD-STM Analysis Complete!\n');
fprintf('========================================\n');
fprintf('Results saved in: %s\n', pwd);
fprintf('Log file: %s_LOGfile.txt\n\n', LOGfile);


%% =========================================================================
%% HELPER FUNCTIONS
%  =========================================================================

function out = iif(condition, true_val, false_val)
    % Inline if function for cleaner conditional expressions
    if condition
        out = true_val;
    else
        out = false_val;
    end
end

