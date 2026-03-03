%  TRUNK SCRIPT: Real Data MT-SBD-STM Analysis
%  ========================================================================
%  Multi-kernel Tensor Shifted Blind Deconvolution for Scanning Tunneling 
%  Microscopy real-data loading, preprocessing, decomposition, and
%  visualization.
%
%  This script follows the UBC LAIR standardized block structure for
%  reproducibility and comprehensive logging. It is the standardized
%  entry point for real-data analysis (the legacy script
%  `MTSBD_block_realdata1.m` remains available for reference).
%
%  WORKFLOW OVERVIEW (CHECKPOINT STAGES):
%  ======================================
%  raw       : load .3ds, copy raw into project, basic metadata
%  preprocess: denoise, crop, slice selection, streak removal, masking, etc.
%  pre-run   : reference slice selection, kernel initialization, ref-slice MT-SBD
%  run       : all-slice MT-SBD (block run)
%  post-run  : visualization and derived analyses
%
%  Each major phase is handled by a wrapper function in `Dong_func/wrapper`
%  that is responsible for logging and checkpointing.
%
%  NOTE: This script is **interactive by default**. Non-interactive (config-
%  only) execution is intended for exact reproduction workflows.
%
% =========================================================================
%% SECTION 0: Path Initialization and Config
% =========================================================================

clc; clear; close all;

% Add paths (run init_sbd from parent directory, mirroring run_synthetic_data.m)
if exist('../init_sbd.m', 'file')
    run('../init_sbd.m');
else
    error('init_sbd.m not found. Please run from scripts/ directory.');
end

% Initialize or preserve core structs
if ~exist('log', 'var');    log = struct();    end
if ~exist('data', 'var');   data = struct();   end
if ~exist('params', 'var'); params = struct(); end
if ~exist('meta', 'var');   meta = struct();   end

% Initialize/upgrade config schema (checkpoint-first design)
% In the future, cfg will primarily come from the loaded parent node.
cfg = init_config();


%% =========================================================================
%% LR01A: Load-Real-01-A; Load .3ds data and create project
%  =========================================================================
%  Loads a real STM/QPI dataset from a .3ds file, creates/opens a project
%  for this dataset, copies the raw data into the project, and writes a
%  "raw" checkpoint (stage = raw).
%
%  Dependencies: load3dsall.m, createProjectStructure (or real-data variant),
%                logUsedBlocks.m, future wrapper loadRealDataset.m
%
% -------------------------------------------------------------------------
% PRESETS: User-configurable parameters
% -------------------------------------------------------------------------
cfg.load.data_file       = 'QPImap012.3ds';  % raw file name (or [] to pick via UI)
cfg.load.smoothing_sigma = 10;               % smoothing for current data

%%%%%%%%%%%%%%%%%% DO NOT EDIT BELOW %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
[log, data, params, meta, cfg] = loadRealDataset(log, data, params, meta, cfg);


%% =========================================================================
%% PR01A: Preprocess-Real-01-A; Preprocess real data
%  =========================================================================
%  Applies preprocessing to the loaded real-data volume: background
%  normalization, Bragg removal, cropping, slice selection, local streak
%  removal (manual/auto), defect masking, streak healing, and final
%  normalization Y used for MT-SBD.
%
%  Dependencies: normalizeBackgroundToZeroMean3D.m, removeBragg.m,
%                maskSquare.m, gridCropMask.m, streak_correction.m,
%                removeLocalStreaks.m, interpolateLocalStreaks.m,
%                gaussianMaskDefects.m, thresholdDefects.m,
%                RemoveStreaks.m, heal_streaks.m,
%                future wrapper preprocessRealData.m
%
% -------------------------------------------------------------------------
% PRESETS: User-configurable parameters
% -------------------------------------------------------------------------
params.preprocessing.interactive     = true;   % use interactive UIs by default
params.preprocessing.save_checkpoint = true;   % write preprocess checkpoint

%%%%%%%%%%%%%%%%%% DO NOT EDIT BELOW %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[log, data, params, meta, cfg] = preprocessRealData(log, data, params, meta, cfg);


%% =========================================================================
%% RS01A: RefSlice-Real-01-A; Decompose reference slice
%  =========================================================================
%  Selects a reference slice, initializes kernels on that slice, and runs
%  MT-SBD to obtain optimal activations and refined kernels. Writes a
%  "pre-run" checkpoint capturing reference-slice results.
%
%  Dependencies: initialize_kernels.m, MT_SBD.m, visualizeRealResult.m,
%                enforce_kernel_polarity.m,
%                future wrapper decomposeRefSliceReal.m
%
% -------------------------------------------------------------------------
% PRESETS: User-configurable parameters
% -------------------------------------------------------------------------
params.refSlice.interactive     = true;   % interactive reference-slice & kernel selection
params.refSlice.save_checkpoint = true;

cfg.reference.same_size       = true;
cfg.reference.kerneltype      = 'selected';
cfg.reference.window_type     = {'gaussian', 2.5};
cfg.reference.square_size     = [80, 80];
cfg.reference.default_ref_slice    = [];   % or set to use as default
cfg.reference.default_num_kernels = [];   % or set to use as default

%%%%%%%%%%%%%%%%%% DO NOT EDIT BELOW %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
[log, data, params, meta, cfg] = decomposeRefSliceReal(log, data, params, meta, cfg);


%% =========================================================================
%% IP01R: Initialize-Proliferation-01-R; Initialize kernels for all slices
%  =========================================================================
%  Uses reference-slice kernel centers/sizes to initialize kernels across
%  all energy slices via proliferation. Writes or updates a "pre-run"
%  checkpoint containing A1_all / A1_all_matrix.
%
%  Dependencies: initialize_kernels_proliferation.m, enforce_kernel_polarity.m,
%                future wrapper initProliferationReal.m
%
% -------------------------------------------------------------------------
% PRESETS: User-configurable parameters
% -------------------------------------------------------------------------
params.proliferation.interactive_size_adjust = false;
params.proliferation.use_matrix_format       = true;

cfg.blockInit.use_matrix  = true;
cfg.blockInit.change_size = false;

%%%%%%%%%%%%%%%%%% DO NOT EDIT BELOW %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
[log, data, params, meta, cfg] = initProliferationReal(log, data, params, meta, cfg);


%% =========================================================================
%% BR01A: Block-Run-Real-01-A; Decompose all slices (block run)
%  =========================================================================
%  Runs MT-SBD on all slices simultaneously (or via an X-regulated variant),
%  using trusted-slice weights and reference-slice initializations. Writes a
%  "run" checkpoint with kernels, activations, bias terms, and extras.
%
%  Dependencies: MTSBD_all_slice.m, MTSBD_Xregulated_all_slices.m,
%                build_auto_trusted_slice_weights.m,
%                future wrapper runAllSlicesReal.m
%
% -------------------------------------------------------------------------
% PRESETS: User-configurable parameters
% -------------------------------------------------------------------------
params.blockRun.interactive     = true;
params.blockRun.save_checkpoint = true;

cfg.blockRun.trusted_ratio_threshold_default = 1.5;
cfg.blockRun.use_default_manual_trusted_slices = true;
cfg.blockRun.show_trusted_plot  = true;
cfg.blockRun.miniloop_iteration = 1;
cfg.blockRun.outerloop_maxIT    = 8;
cfg.blockRun.lambda1_base       = [0.018, 0.016, 0.02, 0.02, 0.015];
cfg.blockRun.phase2             = false;
cfg.blockRun.kplus_factor        = 0.5;
cfg.blockRun.lambda2            = [2e-2, 2e-2];
cfg.blockRun.nrefine            = 3;
cfg.blockRun.signflip           = 0.2;
cfg.blockRun.xpos               = true;
cfg.blockRun.getbias            = true;
cfg.blockRun.Xsolve              = 'FISTA';
cfg.blockRun.use_Xregulated     = false;
cfg.blockRun.allow_custom_update_order = true;

%%%%%%%%%%%%%%%%%% DO NOT EDIT BELOW %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
[log, data, params, meta, cfg] = runAllSlicesReal(log, data, params, meta, cfg);


%% =========================================================================
%% VR01R: Visualize-Real-01-R; Visualize results
%  =========================================================================
%  Optional visualization block for kernels, activations, reconstructions,
%  QPI, and quality metrics. Typically does not create new checkpoints, but
%  may log visualization actions.
%
%  Dependencies: visualizeRealResult.m, d3gridDisplay.m, qpiCalculate.m,
%                future wrapper visualizeRealRun.m
%
% -------------------------------------------------------------------------
% PRESETS: User-configurable parameters
% -------------------------------------------------------------------------
params.visualize.all_slices   = true;
params.visualize.save_figures = false;

%%%%%%%%%%%%%%%%%% DO NOT EDIT BELOW %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
[log, data, params, meta, cfg] = visualizeRealRun(log, data, params, meta, cfg);


%% =========================================================================
%% FINAL: Real-data analysis complete
%  =========================================================================
fprintf('========================================\n');
fprintf('Real-data MT-SBD-STM analysis script finished.\n');
fprintf('========================================\n');

