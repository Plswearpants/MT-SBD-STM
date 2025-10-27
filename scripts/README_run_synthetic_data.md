# MT-SBD-STM Synthetic Data Script - Reorganization Documentation

## Overview

This document describes the reorganized synthetic data script (`run_synthetic_data.m`) that follows the UBC LAIR standardization framework for reproducibility and comprehensive logging.

## Script Structure

### Block-Based Architecture

The script is organized into 8 main functional blocks, each with a unique 5-character identifier following the `ABXXZ` convention:

| Block ID | Name | Purpose | Key Functions |
|----------|------|---------|---------------|
| **GD01A** | Generate-Data-01-A | Generate synthetic 3D STM observation | `properGen_full`, `normalizeBackgroundToZeroMean3D`, `proj2oblique` |
| **IN01A** | Initialize-kNernel-01-A | Initialize kernels for reference slice | `initialize_kernels` |
| **DS01A** | Decompose-Slice-01-A | Decompose reference slice using MT-SBD | `MTSBD_synthetic`, `Asolve_Manopt_tunable`, `Xsolve_FISTA_tunable` |
| **IS01A** | Isolation-Selection-01-A | Find most isolated defect points | `alignActivationMaps`, `padKernels` |
| **IP01A** | Initialize-Proliferation-01-A | Initialize kernels for all slices | `initialize_kernels_proliferation` |
| **DA01A** | Decompose-All-01-A | Decompose all slices simultaneously | `MTSBD_synthetic_all_slice`, `convfft3` |
| **VR01A** | Visualize-Results-01-A | Comprehensive result visualization | `visualizeResults` |
| **WS01A** | Write-Save-01-A | Save workspace and results | Standard MATLAB save |

## Block Identifier Convention

Each block follows the `ABXXZ` format:

- **A (Category)**: G=Generate, I=Initialize, D=Decompose, V=Visualize, W=Write, S=Select, P=Process
- **B (Subcategory)**: D=Data, N=kNernel, S=Slice, P=Proliferation, R=Results, W=Workspace, A=Analysis
- **XX (Number)**: 01, 02, 03, ... (sequential within category)
- **Z (Variant)**: A, B, C, ... (variations of same functionality)

## Standardized Block Structure

Each block follows this pattern:

```matlab
%% =========================================================================
%% ABXXZ: Category-Subcategory-XX-Z; Block Description
%  =========================================================================
%  Edited by [Author], [Date]
%
%  [Detailed description]
%  Dependencies: [list of functions]

% -------------------------------------------------------------------------
% PRESETS: User-configurable parameters
% -------------------------------------------------------------------------
parameter1 = value1;                % Description
parameter2 = value2;                % Description
optional_param = [];                % Optional parameter

%%%%%%%%%%%%%%%%%% DO NOT EDIT BELOW %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% LOG: Block start
LOGcomment = sprintf("Parameter summary...");
LOGcomment = logUsedBlocks(LOGpath, LOGfile, "ABXXZ", LOGcomment, 0);

% [Block execution code]

% LOG: Additional details
LOGcomment = sprintf("Execution details...");
LOGcomment = logUsedBlocks(LOGpath, LOGfile, "  ^  ", LOGcomment, 0);

% Clear preset variables
clearvars parameter1 parameter2 optional_param
```

## Logging System

### Log File Structure

Every execution creates a timestamped log file: `synthetic_run_YYYYMMDD_HHMMSS_LOGfile.txt`

### Logging Pattern

1. **Block Entry**: Block ID + parameter summary
2. **Continuation**: `"  ^  "` + execution details
3. **Results**: `"  ^  "` + key metrics/outputs

### Example Log Entry

```
DATE and TIME         BLOCK   COMMENT
2025-10-27 14:32:15  INIT    === Synthetic Data Analysis Started ===
2025-10-27 14:32:16  GD01A   SNR=5, N_obs=50, resolution=3, density=0.01, slices=2
2025-10-27 14:32:45    ^     Generated Y: 150x150x2, Kernels: 3, Time: 28.53s
2025-10-27 14:32:48    ^     Reference slice: 1
```

## Key Features

### 1. Comprehensive Parameter Documentation

- All user-configurable parameters are in PRESETS section
- Clear comments explain each parameter
- Optional parameters explicitly marked

### 2. Complete Logging

- Every operation logged with timestamp
- Full reproducibility from log file
- Parameters, timing, and results tracked

### 3. Modular Design

- Self-contained blocks
- Clear dependencies listed
- Easy to modify or skip blocks

### 4. Variable Management

- `clearvars` at end of each block prevents namespace pollution
- Structured data storage in `params`, `ref_results`, `allslice_results`, etc.
- No variable carryover between blocks

### 5. User Interaction

- Interactive slice selection with validation
- Progress reporting with timing
- Quality metrics displayed

## Parameter Categories

### Data Generation (GD01A)

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `SNR` | float | 5 | Signal-to-noise ratio |
| `N_obs` | int | 50 | Observation lattice size (pixels) |
| `observation_resolution` | int | 3 | Pixels per lattice site |
| `defect_density` | float | 0.01 | Surface defect density (0-1) |
| `num_slices` | int | 2 | Number of energy slices |
| `LDoS_path` | string | ... | Path to LDoS simulation data |

### Kernel Initialization (IN01A)

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `kernel_selection_type` | string | 'selected' | 'selected' or 'random' |
| `window_type` | cell/string | {'gaussian', 2.5} | Window function type |

### Decomposition (DS01A)

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `use_Xregulated` | bool | false | Use X-regularized version |
| `initial_iteration` | int | 1 | Inner solver iterations (start) |
| `maxIT` | int | 15 | Outer alternating iterations |
| `lambda1` | array | [3e-2, ...] | L1 regularization (Phase I) |
| `phase2_enable` | bool | false | Enable Phase II refinement |
| `lambda2` | array | [1e-2, ...] | L1 regularization (Phase II) |
| `nrefine` | int | 5 | Number of refinement steps |
| `xpos` | bool | true | Enforce positive activations |
| `getbias` | bool | true | Extract constant bias |
| `Xsolve_method` | string | 'FISTA' | 'FISTA' or 'pdNCG' |

### Isolation Analysis (IS01A)

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `isolation_threshold_factor` | float | 10 | Threshold = max/factor |
| `target_kernel_size_type` | string | 'kernel_sizes_all' | Size selection method |

### Proliferation (IP01A)

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `window_type_proliferation` | cell | {'gaussian', 2.5} | Window function |
| `interactive_size_adjust` | bool | false | Allow interactive resizing |
| `use_matrix_format` | bool | true | Matrix vs cell output |

### All-Slice Decomposition (DA01A)

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `use_Xregulated_allslice` | bool | false | Use X-regularized version |
| `use_reference_init` | bool | true | Initialize from reference |
| `maxIT_allslice` | int | 15 | Max iterations |

### Visualization (VR01A)

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `visualize_all_slices` | bool | true | Visualize all vs ref only |
| `save_figures` | bool | false | Save figure files |

### Workspace Save (WS01A)

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `save_workspace` | bool | true | Save complete workspace |
| `save_results_only` | bool | false | Save key results only |
| `results_path` | string | pwd | Save directory |

## Function Dependencies

### By Block

**GD01A (Generate Data)**
- `properGen_full.m` - Main synthetic data generator
- `normalizeBackgroundToZeroMean3D.m` - Normalization
- `proj2oblique.m` - Manifold projection
- `d3gridDisplay.m` - 3D visualization

**IN01A (Initialize Kernels)**
- `initialize_kernels.m` - Interactive kernel selection
- `apply_window.m` - Window function application

**DS01A (Decompose Slice)**
- `MTSBD_synthetic.m` - Main SBD algorithm
- `MTSBD_synthetic_Xregulated.m` - X-regularized variant
- `Asolve_Manopt_tunable.m` - Kernel solver
- `Xsolve_FISTA_tunable.m` - Activation solver
- `computeQualityMetrics.m` - Quality assessment
- `showims.m` - Progress visualization
- `update_config.m` - Solver configuration

**IS01A (Isolation Selection)**
- `alignActivationMaps.m` - Activation alignment
- `padKernels.m` - Kernel padding

**IP01A (Initialize Proliferation)**
- `initialize_kernels_proliferation.m` - Proliferation initialization

**DA01A (Decompose All)**
- `MTSBD_synthetic_all_slice.m` - Multi-slice SBD
- `MTSBD_synthetic_Xregulated_all_slices.m` - X-regularized variant
- `convfft3.m` - 3D convolution

**VR01A (Visualize Results)**
- `visualizeResults.m` - Comprehensive visualization

## Comparison with Original Script

| Aspect | Original (`MTSBD_block_synthetic.m`) | New (`run_synthetic_data.m`) |
|--------|--------------------------------------|------------------------------|
| **Structure** | 3 informal blocks with comments | 8 formal blocks with IDs |
| **Logging** | None | Comprehensive with timestamps |
| **Parameters** | Scattered throughout | Organized in PRESETS sections |
| **Documentation** | Minimal inline comments | Detailed headers and descriptions |
| **Reproducibility** | Manual tracking required | Automatic via log files |
| **Variable Management** | Variables persist | Cleaned after each block |
| **Modularity** | Monolithic flow | Independent, reusable blocks |
| **Error Handling** | Basic validation | Comprehensive validation |

## Usage Example

### Basic Usage

```matlab
% Navigate to scripts directory
cd scripts/

% Run the script
run_synthetic_data

% Follow prompts:
% 1. Wait for data generation
% 2. View 3D grid display and enter reference slice number
% 3. Select kernel regions interactively
% 4. Script runs automatically through decomposition
% 5. Results are visualized and saved
```

### Customizing Parameters

```matlab
% Edit the PRESETS sections in each block before running
% For example, in Block GD01A:
SNR = 10;                    % Increase SNR
N_obs = 100;                 % Larger observation
defect_density = 0.005;      % Lower density

% Then run as normal
```

### Skipping Blocks

To skip a block (e.g., if data already generated):

1. Comment out the entire block from `%%` to `clearvars`
2. Load required variables from previous run
3. Continue with subsequent blocks

## Future Enhancements

### Planned Improvements

1. **Configuration Files**: Move parameters to external config files
2. **Batch Processing**: Add support for parameter sweeps
3. **GUI**: Optional graphical interface for parameter selection
4. **Resume Capability**: Checkpoint system for long runs
5. **Parallel Processing**: Multi-core support for proliferation
6. **Performance Profiling**: Built-in timing and memory tracking

### Additional Blocks (Future)

- **PM01A**: Process-Metrics-01-A - Advanced metric computation
- **ST01A**: Select-Threshold-01-A - Adaptive threshold selection
- **GP01A**: Generate-Plots-01-A - Publication-quality figures
- **WR01A**: Write-Report-01-A - Automated report generation

## Troubleshooting

### Common Issues

**Issue**: Log file not created
- **Solution**: Check write permissions in scripts/ directory
- **Solution**: Verify `logUsedBlocks.m` is in path

**Issue**: Function not found
- **Solution**: Run `init_sbd.m` from parent directory first
- **Solution**: Check that all dependencies are in correct folders

**Issue**: Reference slice selection fails
- **Solution**: Enter number between 1 and num_slices
- **Solution**: Ensure 3D display is visible

**Issue**: Out of memory during all-slice decomposition
- **Solution**: Reduce `N_obs` or `num_slices`
- **Solution**: Use `save_results_only = true` in WS01A

## References

- UBC LAIR Template: `template/TEMPLATE_Script_DOCUMENTATION.md`
- Quick Reference: `template/TEMPLATE_script_QUICK_REFERENCE.md`
- Logging Function: `Dong_func/basic/logUsedBlocks.m`
- Original Script: `scripts/MTSBD_block_synthetic.m`

## Version History

- **v1.0** (2025-10-27): Initial standardized implementation
  - 8 formal blocks with comprehensive logging
  - Full parameter documentation
  - UBC LAIR compliance

---

**Author**: [Your Name]  
**Last Updated**: 2025-10-27  
**Project**: MT-SBD-STM

