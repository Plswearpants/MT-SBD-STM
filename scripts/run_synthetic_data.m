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
%  WORKFLOW OVERVIEW (THREE-PHASE STRUCTURE):
%  ===========================================
%  PHASE 1: PRE-RUN (Synthetic Data Generation & Initialization)
%    Block 1 (GD01A): Generate synthetic test data + auto kernel initialization
%    Block 2 (IK01A): [OPTIONAL] Manual kernel initialization
%    → Auto-save dataset (auto.mat or manualXX.mat)
%
%  PHASE 2: POST-RUN (Algorithm Execution)
%    Block 3 (DS01A): Decompose reference slice (find optimal activation)
%    Block 4 (IS01A): Find most isolated points for kernel proliferation
%    Block 5 (IP01A): Initialize kernels for all slices (proliferation)
%    Block 6 (DA01A): Decompose all slices simultaneously
%    → Auto-save run results (runXX/runXX.mat)
%
%  PHASE 3: VISUALIZATION
%    Block 7 (VR01A): Visualize results
%    (No auto-save for visualization)
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
%% SECTION 0: Path Initialization
% =========================================================================

% Clear workspace (optional - comment out if you want to keep variables)
clc; clear; close all;

% Add paths (run init_sbd from parent directory)
if exist('../init_sbd.m', 'file')
    run('../init_sbd.m');
else
    error('init_sbd.m not found. Please run from scripts/ directory.');
end

% Note: Logging is initialized in GD01A (dataset-specific log file)

%% =========================================================================
%% GD01A: Generate-Data-01-A; Generate Normalized Synthetic Test Data
%  =========================================================================
%  Edited by Dong Chen, 2025-10-27
%
%  Generates synthetic 3D STM observation data with multiple kernels across
%  energy slices. Uses LDoS simulation data and creates realistic defect
%  patterns with controlled SNR.
%
%  Dependencies: generateSyntheticData.m

% -------------------------------------------------------------------------
% PRESETS: User-configurable parameters
% -------------------------------------------------------------------------
params.synGen.SNR = 5;                        % Signal-to-noise ratio
params.synGen.N_obs = 50;                     % Observation lattice size (pixels)
params.synGen.observation_resolution = 3;     % Resolution: pixels per lattice site
params.synGen.defect_density = 0.01;          % Surface defect density (0-1)
params.synGen.num_slices = 2;                 % Number of energy slices
params.synGen.LDoS_path = 'C:\Users\CAD\Documents\GitHub\MT-SBD-STM\examples\example_data\LDoS_single_defects_self=0.6_save.mat'; % make this an UI option 
params.synGen.vis_generation = false;         % Show intermediate generation steps
params.synGen.normalization_type = 'dynamic'; % 'dynamic' or 'static'
params.synGen.ref_slice = [];                 % Reference slice ([] for interactive selection)

%%%%%%%%%%%%%%%%%% DO NOT EDIT BELOW %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create project structure first (fresh project for each data generation)
fprintf('Creating project structure...\n');
meta = createProjectStructure();
fprintf('Project folder: %s\n\n', meta.project_name);

% Create fresh log file for this dataset (prevents log contamination)
% Log file will be saved with the dataset
fprintf('Initializing dataset log file...\n');
dataset_log_name = sprintf('dataset_%s', meta.timestamp);
log.path = meta.project_path;
log.file = dataset_log_name;
LOGcomment = sprintf("Dataset generation session: %s", meta.project_name);
LOGcomment = logUsedBlocks(log.path, log.file, "GD01A", LOGcomment, 1);
fprintf('Dataset log file: %s_LOGfile.txt\n\n', log.file);

% Generate synthetic data (with internal logging)
fprintf('Generating synthetic test data...\n');
[data, params] = generateSyntheticData(log, params);
fprintf('Data generation complete.\n\n');

% Auto-initialize kernels from ground truth (with internal logging)
fprintf('Auto-initializing kernels from ground truth...\n');
[data, params] = autoInitializeKernels(log, data, params);
fprintf('Auto initialization complete.\n\n');

% Auto-save dataset with auto initialization
fprintf('Saving dataset with auto initialization...\n');
params = organizeParams(params, 'write');  % Convert to hierarchical for storage
meta = saveDataset(log, data, params, meta);


%% =========================================================================
%% IK01A: [OPTIONAL] Initialize-Kernel-01-A; Manual Kernel Initialization
%  =========================================================================
%  Edited by Dong Chen, 2025-11-04
%
%  OPTIONAL BLOCK: Comment out this entire block to skip manual initialization.
%  If run, this will overwrite the auto initialization and save as manualXX.mat
%
%  Interactive manual kernel initialization from the reference slice.
%  User selects regions containing isolated defect patterns to use as
%  initial kernel guesses.
%
%  Dependencies: initializeKernelsRef.m
%
%  TO SKIP THIS BLOCK: Select from %% IK01A to the next %% and comment out (Ctrl+R)

% -------------------------------------------------------------------------
% PRESETS: User-configurable parameters
% -------------------------------------------------------------------------
params.initialization.kernel_selection_type = 'selected';     % 'selected' or 'random'
params.initialization.window_type = {'gaussian', 2.5};        % Window function for kernels
                                                        % Options: 'hann', 'hamming', 'blackman'
                                                        %          {'gaussian', alpha}
                                                        %          {'kaiser', beta}
                                                        %          {} for no window

%%%%%%%%%%%%%%%%%% DO NOT EDIT BELOW %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize kernels manually (overwrites auto initialization, with internal logging)
fprintf('Running manual kernel initialization...\n');
[data, params] = initializeKernelsRef(log, data, params);

% Auto-save dataset with manual initialization
fprintf('Saving dataset with manual initialization...\n');
params = organizeParams(params, 'write');  % Convert to hierarchical for storage
meta = saveDataset(log, data, params, meta);
fprintf('Dataset saved: %s\n', meta.dataset_file);
fprintf('Location: %s\n\n', meta.project_path);


%% =========================================================================
%% LD01A: Load-Dataset-01-A; Load Dataset for Algorithm Run
%  =========================================================================
%  Edited by Dong Chen, 2025-11-03
%
%  Loads a pre-generated dataset (with auto or manual initialization) to
%  prepare for algorithm execution. This block opens a file browser UI where
%  the user selects a .mat file. The function automatically finds and verifies
%  the corresponding log file in the same directory.
%
%  Prerequisites: Run GD01A to create project structure
%  Dependencies: loadWorkspace.m
%
%  NOTE: This block has no user-configurable presets (interactive selection).

%%%%%%%%%%%%%%%%%% DO NOT EDIT BELOW %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load dataset interactively (file browser UI) with memory protection
% User selects .mat file, function automatically finds corresponding log file
% Handle case where meta might not exist (e.g., jumping directly to this block)
if exist('meta', 'var') && isstruct(meta)
    current_meta = meta;
else
    current_meta = struct();
end

if exist('log', 'var') && exist('data', 'var') && exist('params', 'var')
    [log, data, params, meta] = loadWorkspace('current_log', log, 'current_data', data, 'current_params', params, 'current_meta', current_meta, 'title', 'LOAD DATASET FOR ALGORITHM RUN', 'track_selection', true);
else
    [log, data, params, meta] = loadWorkspace('current_meta', current_meta, 'title', 'LOAD DATASET FOR ALGORITHM RUN', 'track_selection', true);
end

%% =========================================================================
%% DS01A: Decompose-Slice-01-A; Decompose Reference Slice
%  =========================================================================
%  Edited by Dong Chen, 2025-10-28
%
%  Runs the MT-SBD algorithm on the reference slice to find optimal
%  activation maps and refine kernel estimates. This provides the starting
%  point for multi-slice analysis.
%
%  Prerequisites: Run LD01A to load a dataset first
%  Dependencies: decomposeReferenceSlice.m

% -------------------------------------------------------------------------
% PRESETS: User-configurable parameters
% -------------------------------------------------------------------------
% Phase I settings
params.mcsbd_slice.initial_iteration = 1;          % Manopt/FISTA inner iterations (start)
params.mcsbd_slice.maxIT = 15;                     % Number of outer alternating iterations
params.mcsbd_slice.lambda1 = 3e-2;                 % L1 regularization (Phase I) - scalar or vector

% Phase II settings (refinement)
params.mcsbd_slice.phase2_enable = false;          % Enable Phase II refinement
params.mcsbd_slice.lambda2 = 1e-2;                 % L1 regularization (Phase II final) - scalar or vector
params.mcsbd_slice.nrefine = 5;                    % Number of refinement steps
params.mcsbd_slice.kplus_factor = 0.5;             % Sphere lifting padding factor

% Algorithm parameters
params.mcsbd_slice.signflip_threshold = 0.2;       % Sign flip detection threshold
params.mcsbd_slice.xpos = true;                    % Enforce positive activations
params.mcsbd_slice.getbias = true;                 % Extract constant bias term
params.mcsbd_slice.Xsolve_method = 'FISTA';        % 'FISTA' or 'pdNCG'

% Initialization options
params.mcsbd_slice.use_xinit = [];                 % Initial X guess ([] for none)

% Display options
params.mcsbd_slice.show_progress = true;           % Show progress during optimization

%%%%%%%%%%%%%%%%%% DO NOT EDIT BELOW %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Decompose reference slice (wrapper handles logging internally)
fprintf('Decomposing reference slice...\n');
[data, params] = decomposeReferenceSlice(log, data, params);

% Clear preset variables


%% =========================================================================
%% IS01A: Isolation-Selection-01-A; Find Most Isolated Points
%  =========================================================================
%  Edited by Dong Chen, 2025-10-29
%
%  Analyzes the activation maps from the reference slice to find the most
%  isolated defect for each kernel. These points are used as centers for
%  kernel initialization across all energy slices.
%
%  Dependencies: findIsolatedPoints.m

% -------------------------------------------------------------------------
% PRESETS: User-configurable parameters
% -------------------------------------------------------------------------
isolation_threshold_factor = 10;    % Threshold = max/factor for defect detection
target_kernel_size_type = 'kernel_sizes_all';  % 'ref_kernel_sizes',
                                               % 'kernel_sizes_cap',
                                               % 'kernel_sizes_all'
show_distributions = true;          % Show activation value distributions
show_isolation_maps = true;         % Show isolation analysis visualizations

%%%%%%%%%%%%%%%%%% DO NOT EDIT BELOW %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Find most isolated points for 3D initialization (wrapper handles logging internally)
fprintf('Finding most isolated points for 3D initialization...\n');
[data, params, isolation_results] = findIsolatedPoints(log, data, params, ...
    'isolation_threshold_factor', isolation_threshold_factor, ...
    'target_kernel_size_type', target_kernel_size_type, ...
    'show_distributions', show_distributions, ...
    'show_isolation_maps', show_isolation_maps);

fprintf('Isolation selection complete.\n\n');

% Clear preset variables
clearvars isolation_threshold_factor target_kernel_size_type
clearvars show_distributions show_isolation_maps


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
LOGcomment = logUsedBlocks(log.path, log.file, "IP01A", LOGcomment, 0);

fprintf('Initializing kernels for all %d slices using proliferation...\n', params.num_slices);

% Convert most_isolated_points to matrix format
kernel_centers = zeros(params.num_kernels, 2);
for k = 1:params.num_kernels
    if ~isempty(data.mcsbd_slice.most_isolated_points{k})
        kernel_centers(k,:) = data.mcsbd_slice.most_isolated_points{k};
    else
        error('No isolated point found for kernel %d', k);
    end
end

% Initialize kernels for all slices
A1_all = cell(params.num_slices, params.num_kernels);

for s = 1:params.num_slices
    fprintf('  Slice %d/%d...', s, params.num_slices);
    
    if strcmp(params.target_kernel_size_type, 'kernel_sizes_all') || ndims(params.target_kernel_sizes) == 3
        % Use slice-specific kernel sizes
        target_sizes_slice = squeeze(params.kernel_sizes(s,:,:));
    else
        % Use uniform kernel sizes
        target_sizes_slice = params.target_kernel_sizes;
    end
    
    A1_all(s,:) = initialize_kernels_proliferation(data.synGen.Y(:,:,s), params.num_kernels, ...
        kernel_centers, window_type_proliferation, target_sizes_slice, ...
        'interactive', interactive_size_adjust);
    
    fprintf(' done.\n');
end

% Visualize initialized kernels for reference slice
figure('Name', 'IP01A: Initialized Kernels (Reference Slice)');
for k = 1:params.num_kernels
    subplot(1, params.num_kernels, k);
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
init_results.kernel_sizes = params.target_kernel_sizes;

fprintf('Kernel initialization complete for all slices.\n\n');

% Convert to matrix format if requested
if use_matrix_format
    A1_all_matrix = cell(params.num_kernels, 1);
    for k = 1:params.num_kernels
        A1_all_matrix{k} = zeros(size(A1_all{1,k},1), size(A1_all{1,k},2), params.num_slices);
        for s = 1:params.num_slices
            A1_all_matrix{k}(:,:,s) = A1_all{s,k};
        end
    end
    fprintf('Converted to matrix format: %d kernels x [%dx%dx%d]\n\n', ...
        params.num_kernels, size(A1_all_matrix{1},1), size(A1_all_matrix{1},2), params.num_slices);
end

% LOG: Initialization complete
LOGcomment = sprintf("Initialized %d kernels for %d slices at centers: %s", ...
    params.num_kernels, params.num_slices, mat2str(kernel_centers));
LOGcomment = logUsedBlocks(log.path, log.file, "  ^  ", LOGcomment, 0);

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
% Use reference slice results as initialization
use_reference_init = true;          % Initialize X with reference results

% Display options
show_allslice_progress = true;      % Show progress during optimization

% Algorithm parameters inherited from DS01A
maxIT_allslice = 15;                % Max iterations for all-slice decomposition

%%%%%%%%%%%%%%%%%% DO NOT EDIT BELOW %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% LOG: Block start
LOGcomment = sprintf("Init from ref: %d, maxIT=%d", use_reference_init, maxIT_allslice);
LOGcomment = logUsedBlocks(log.path, log.file, "DA01A", LOGcomment, 0);

fprintf('Running MT-SBD decomposition on all %d slices...\n', params.num_slices);

% Prepare kernel sizes and observation
kernel_sizes_used = squeeze(max(params.kernel_sizes, [], 1));
Y_used = data.synGen.Y;
A1_used = A1_all_matrix;

% Set up display functions
if show_allslice_progress
    figure('Name', 'DA01A: All-Slice Decomposition Progress');
    dispfun_allslice = cell(1, params.num_kernels);
    for n = 1:params.num_kernels
        dispfun_allslice{n} = @(Y, A, X, kernel_sizes, kplus) ...
            showims(Y_used, A1_used{n}, data.synGen.X0_ref(:,:,n), A, X, kernel_sizes, kplus, 1);
    end
else
    dispfun_allslice = cell(1, params.num_kernels);
    for n = 1:params.num_kernels
        dispfun_allslice{n} = @(Y, A, X, kernel_sizes, kplus) 0;
    end
end

% Update parameters for all-slice processing
params_allslice = ref_results.params;
params_allslice.X0 = data.synGen.X0;
params_allslice.A0 = data.mcsbd_slice.A0_used;

% Initialize X from reference slice if requested
if use_reference_init
    params_allslice.xinit = cell(1, params.num_kernels);
    for k = 1:params.num_kernels
        params_allslice.xinit{k}.X = data.mcsbd_slice.X(:,:,k);
        b_temp = ref_results.extras.phase1.biter(k);
        params_allslice.xinit{k}.b = repmat(b_temp, [params.num_slices, 1]);
    end
else
    params_allslice.xinit = [];
end

% Run all-slice decomposition
tic;
[Aout_slice, Xout_slice, bout_slice, slice_extras] = MTSBD_synthetic_all_slice(...
    Y_used, kernel_sizes_used, params_allslice, dispfun_allslice, A1_used, initial_iteration, maxIT_allslice);
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
LOGcomment = logUsedBlocks(log.path, log.file, "  ^  ", LOGcomment, 0);

% Clear preset variables
clearvars use_reference_init show_allslice_progress maxIT_allslice
clearvars kernel_sizes_used Y_used A1_used dispfun_allslice
clearvars params_allslice b_temp allslice_time n k


%% =========================================================================
%% PHASE 2 COMPLETION: Save Run Results (Post-Run Phase)
%  =========================================================================
%  Auto-save algorithm execution results.
%  Saved in: runXX/runXX.mat

fprintf('\n========================================\n');
fprintf('PHASE 2 COMPLETE: Post-Run\n');
fprintf('========================================\n');

% Auto-save run results
params = organizeParams(params, 'write');  % Convert to hierarchical for storage
meta = saveRun(log, data, params, meta);

fprintf('Run saved: %s\n', meta.run_file);
fprintf('Ready for Phase 3 (Visualization)\n');
fprintf('========================================\n\n');


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
LOGcomment = logUsedBlocks(log.path, log.file, "VR01A", LOGcomment, 0);

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
    fprintf('  Visualizing slice %d/%d...\n', s, params.num_slices);
    
    % Prepare slice-specific data
    Y_slice = data.synGen.Y(:,:,s);
    A0_slice = cell(1, params.num_kernels);
    for k = 1:params.num_kernels
        if iscell(data.mcsbd_slice.A0_used)
            A0_slice{k} = data.mcsbd_slice.A0_used{k}(:,:,s);
        else
            A0_slice{k} = data.mcsbd_slice.A0_used{s,k};
        end
    end
    Aout_slice_cell = cell(1, params.num_kernels);
    for k = 1:params.num_kernels
        Aout_slice_cell{k} = Aout_cell{s,k};
    end
    
    % Call visualizeResults for this slice
    visualizeResults(Y_slice, A0_slice, Aout_slice_cell, data.synGen.X0, Xout_slice, ...
        bout_slice(s,:), slice_extras, [s, 1]);
    
    % Save figures if requested
    if save_figures
        fig_name = sprintf('%s_slice%02d_results', log.file, s);
        savefig(gcf, fullfile(log.path, [fig_name '.fig']));
        LOGcomment = sprintf("Saved figure: %s.fig", fig_name);
        LOGcomment = logUsedBlocks(log.path, log.file, "  ^  ", LOGcomment, 0);
    end
end

fprintf('Visualization complete.\n\n');

% LOG: Visualization complete
LOGcomment = sprintf("Visualized %d slices", length(slices_to_viz));
LOGcomment = logUsedBlocks(log.path, log.file, "  ^  ", LOGcomment, 0);

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
LOGcomment = logUsedBlocks(log.path, log.file, "WS01A", LOGcomment, 0);

if save_workspace || save_results_only
    fprintf('Saving results...\n');
    
    % Create results filename
    results_filename = sprintf('%s_results.mat', log.file);
    results_fullpath = fullfile(results_path, results_filename);
    
    if save_results_only
        % Save only essential results
        save(results_fullpath, 'Y', 'A0', 'X0', 'params', ...
            'ref_results', 'allslice_results', 'isolation_analysis', ...
            'num_kernels', 'num_slices', 'log.file', '-v7.3');
        fprintf('Results saved to: %s\n', results_fullpath);
    else
        % Save complete workspace
        save(results_fullpath, '-v7.3');
        fprintf('Complete workspace saved to: %s\n', results_fullpath);
    end
    
    % Copy log file to results location
    log_source = fullfile(log.path, [log.file '_log.file.txt']);
    log_dest = fullfile(results_path, [log.file '_log.file.txt']);
    if ~strcmp(log_source, log_dest)
        copyfile(log_source, log_dest);
    end
    
    fprintf('Log file saved to: %s\n', log_dest);
    fprintf('\n');
    
    % LOG: Save complete
    LOGcomment = sprintf("Saved to: %s", results_filename);
    LOGcomment = logUsedBlocks(log.path, log.file, "  ^  ", LOGcomment, 0);
end

% Clear preset variables
clearvars save_workspace save_results_only results_path
clearvars results_filename results_fullpath log_source log_dest


%% =========================================================================
%% FINAL: Analysis Complete
%  =========================================================================

% LOG: Final entry
LOGcomment = "=== Synthetic Data Analysis Complete ===";
LOGcomment = logUsedBlocks(log.path, log.file, "DONE ", LOGcomment, 0);

fprintf('========================================\n');
fprintf('MT-SBD-STM Analysis Complete!\n');
fprintf('========================================\n');
fprintf('Results saved in: %s\n', pwd);
fprintf('Log file: %s_log.file.txt\n\n', log.file);


%% =========================================================================
%% HELPER FUNCTIONS
%  =========================================================================
