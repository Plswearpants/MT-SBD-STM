# MT-SBD-STM Parameter Glossary

**Last Updated**: 2025-10-27  
**Author**: Dong Chen  
**Project**: Multi-kernel Tensor Shifted Blind Deconvolution for STM

---

## Table of Contents

1. [Data Generation Parameters](#data-generation-parameters)
2. [Kernel Initialization Parameters](#kernel-initialization-parameters)
3. [Decomposition Parameters](#decomposition-parameters)
4. [Isolation Analysis Parameters](#isolation-analysis-parameters)
5. [Block Initialization Parameters (IB01A)](#block-initialization-parameters-ib01a)
6. [Visualization Parameters](#visualization-parameters)
7. [Data Structures](#data-structures)

---

## Data Generation Parameters

### Block: GD01A (Generate-Data-01-A)
### Function: `generateSyntheticData.m`

| Parameter | Type | Units | Range | Default | Description |
|-----------|------|-------|-------|---------|-------------|
| **SNR** | scalar | - | > 0 | *required* | Signal-to-noise ratio. Higher values = less noise. Typical: 2-10 |
| **N_obs** | integer | pixels | > 0 | *required* | Observation lattice size (square). Determines spatial dimensions of final observation: `N_obs * observation_resolution` |
| **observation_resolution** | integer | pixels/site | > 0 | *required* | Resolution factor: number of pixels per lattice site. Typical: 3-5 |
| **defect_density** | float | - | (0, 1) | *required* | Surface defect density. Fraction of lattice sites with defects. Typical: 0.001-0.05 |
| **num_slices** | integer | - | > 0 | *required* | Number of energy slices in 3D observation. Each slice represents different energy level |
| **LDoS_path** | string | - | - | *required* | Path to LDoS (Local Density of States) simulation data file. Must be .mat file |
| **vis_generation** | logical | - | {true, false} | false | Show intermediate visualization steps during data generation |
| **normalization_type** | string | - | {'dynamic', 'static'} | 'dynamic' | Normalization method for observation. Dynamic adjusts per slice |
| **ref_slice** | integer | - | [1, num_slices] | [] | Reference slice for initial decomposition. Empty = interactive selection |

### Derived/Output Parameters

| Parameter | Type | Units | Description |
|-----------|------|-------|-------------|
| **num_kernels** | integer | - | Number of distinct kernel types detected in LDoS data |
| **kernel_sizes** | array [S×K×2] | pixels | Size of each kernel: [num_slices × num_kernels × (height, width)] |
| **selected_indices** | vector | - | Indices of LDoS slices selected for kernel generation |
| **cutoff_M** | vector | - | Cutoff values used for each kernel during generation |

---

## Kernel Initialization Parameters

### Block: IN01A (Initialize-kNernel-01-A)
### Function: `initializeKernelsRef.m`

| Parameter | Type | Units | Range | Default | Description |
|-----------|------|-------|-------|---------|-------------|
| **kernel_selection_type** | string | - | {'selected', 'random'} | 'selected' | Method for initial kernel selection. 'selected' = interactive, 'random' = random initialization |
| **window_type** | cell/string | - | see below | {'gaussian', 2.5} | Window function applied to kernels to reduce edge effects |

#### Derived/Output Parameters

| Parameter | Type | Units | Description |
|-----------|------|-------|-------------|
| **kernel_sizes_ref** | array [K×2] | pixels | Kernel sizes for reference slice: [num_kernels × (height, width)] |

#### Window Type Options:
- `'hann'` - Hann window (raised cosine)
- `'hamming'` - Hamming window
- `'blackman'` - Blackman window
- `{'gaussian', alpha}` - Gaussian window with parameter alpha (typical: 2-3)
- `{'kaiser', beta}` - Kaiser window with parameter beta (typical: 3-5)
- `{}` or `''` - No window (rectangular)

---

## Decomposition Parameters

### Block: DS01A (Decompose-Slice-01-A)
### Function: `decomposeReferenceSlice.m`

#### Notes
- Lambda parameters can be specified as scalars (applied to all kernels) or vectors (per-kernel values)
- The wrapper auto-sizes lambda vectors to match num_kernels if scalar provided

#### Phase I Settings (Initial Decomposition)

| Parameter | Type | Units | Range | Default | Description |
|-----------|------|-------|-------|---------|-------------|
| **initial_iteration** | integer | iterations | ≥ 1 | 1 | Number of inner solver iterations at start. Increases gradually |
| **maxIT** | integer | iterations | > 0 | 15 | Maximum number of outer alternating iterations |
| **lambda1** | scalar/vector | - | > 0 | 3e-2 | L1 regularization for Phase I. Scalar applies to all kernels, or per-kernel vector |

#### Phase II Settings (Refinement - Optional)

| Parameter | Type | Units | Range | Default | Description |
|-----------|------|-------|-------|---------|-------------|
| **phase2_enable** | logical | - | {true, false} | false | Enable Phase II refinement with sphere lifting |
| **lambda2** | scalar/vector | - | > 0 | 1e-2 | Final L1 regularization for Phase II. Scalar applies to all kernels, or per-kernel vector |
| **nrefine** | integer | steps | ≥ 1 | 5 | Number of refinement steps. Lambda decreases from lambda1 to lambda2 |
| **kplus_factor** | float | - | > 0 | 0.5 | Sphere lifting padding factor: `kplus = kplus_factor × kernel_size` |

#### Algorithm Parameters

| Parameter | Type | Units | Range | Default | Description |
|-----------|------|-------|-------|---------|-------------|
| **signflip_threshold** | float | - | [0, 1] | 0.2 | Threshold for detecting sign flips in optimization |
| **xpos** | logical | - | {true, false} | true | Enforce non-negative constraint on activations X |
| **getbias** | logical | - | {true, false} | true | Extract constant bias term from observation |
| **Xsolve_method** | string | - | {'FISTA', 'pdNCG'} | 'FISTA' | Solver for activation map X. FISTA = Fast Iterative Shrinkage-Thresholding |
| **use_xinit** | struct/[] | - | - | [] | Initial guess for X. Empty = start from scratch |
| **show_progress** | logical | - | {true, false} | true | Display optimization progress during decomposition |

#### Quality Metrics (Output)

| Parameter | Type | Units | Description |
|-----------|------|-------|-------------|
| **activation_metrics** | matrix [IT×K] | - | Activation similarity at each iteration for each kernel |
| **kernel_quality_factors** | matrix [IT×K] | - | Kernel quality score at each iteration for each kernel |

---

## Isolation Analysis Parameters

### Block: IS01A (Isolation-Selection-01-A)
### Function: `findIsolatedPoints.m`

| Parameter | Type | Units | Range | Default | Description |
|-----------|------|-------|-------|---------|-------------|
| **isolation_threshold_factor** | float | - | > 1 | 10 | Defect detection threshold: `threshold = max_activation / factor` |
| **target_kernel_size_type** | string | - | see below | 'kernel_sizes_all' | Method for determining target kernel sizes across slices |
| **show_distributions** | logical | - | {true, false} | true | Display activation value distribution histograms and CDFs |
| **show_isolation_maps** | logical | - | {true, false} | true | Display isolation analysis visualizations |

#### Target Kernel Size Options:
- `'ref_kernel_sizes'` - Use reference slice kernel sizes for all slices
- `'kernel_sizes_cap'` - Use maximum kernel size across all slices
- `'kernel_sizes_all'` - Use slice-specific kernel sizes (maintains energy dependence)

#### Derived/Output Parameters

| Parameter | Type | Units | Description |
|-----------|------|-------|-------------|
| **kernel_centers** | array [K×2] | pixels | Most isolated point for each kernel: [y, x] coordinates |
| **target_kernel_sizes** | array | pixels | Kernel sizes for 3D initialization (depends on size type) |
| **A0_used** | cell | - | Ground truth kernels padded to target sizes if needed |

#### Isolation Results Struct

Contains detailed analysis results:

| Field | Type | Description |
|-------|------|-------------|
| **defect_positions** | cell {1×K} | All detected defect positions [N×2] above threshold per kernel |
| **most_isolated_points** | cell {1×K} | Most isolated point [y, x] per kernel (cell format) |
| **used_most_isolated_points** | cell {1×K} | Isolated points aligned with ground truth |
| **isolation_scores** | cell {1×K} | Distance-based isolation scores for detected defects |
| **num_defects** | vector [1×K] | Number of defects detected per kernel |
| **offset** | array [K×2] | Alignment offset from ground truth activation maps |
| **kernel_centers** | array [K×2] | Centers in matrix format (same as data.kernel_centers) |

---

## Block Initialization Parameters (IB01A)

### Block: IB01A (Initialize-Block-01-A)
### Function: `initializeBlock.m` (unified wrapper); dispatches to `initializeProliferation.m` or block_manual path; core: `initialize_kernels_proliferation.m`

Single entry point for initializing kernels across all slices. **Method** (`block_init_method`) selects how centers are obtained:

| Parameter | Type | Units | Range | Default | Description |
|-----------|------|-------|-------|---------|-------------|
| **block_init_method** | string | - | see below | 'proliferation' | How to obtain block initialization |
| **A1_matrix_unify_size** | string | - | see below | 'max_per_kernel' | How to unify sizes when building A_all_init_matrix (default for block_manual) |
| **window_type_proliferation** | cell | - | see IN01A | {'gaussian', 2.5} | Window function for block-initialized kernels |
| **interactive_size_adjust** | logical | - | {true, false} | false | Allow interactive kernel size adjustment |
| **use_matrix_format** | logical | - | {true, false} | true | Output A_all_init_matrix as {K×1} [H×W×S]; false = cell {S×K} only |
| **show_reference_kernels** | logical | - | {true, false} | true | Display initialized kernels for reference slice |

**block_init_method options:**  
- `'proliferation'` – Use most isolated points from IS01A as centers; initialize kernels at those positions on every slice. Requires IS01A and `target_kernel_sizes`.  
- `'block_manual'` – Always runs manual ref-slice selection (same as IK01A); previous data.slice is not used (manual overwrites it). Centers are derived from the new ref-slice kernels and propagated to all slices. Kernel size default: `A1_matrix_unify_size = 'max_per_kernel'`.

**A1_matrix_unify_size options:**  
- `'max_per_kernel'` – Pad each slice to max height/width over all slices for that kernel (matches DA01A).  
- `'ref_slice'` – Use reference-slice dimensions; pad or center-crop other slices to match.

#### Output (in `data.block`)

| Field | Type | Description |
|-------|------|-------------|
| **A_all_init** | cell {S×K} | Initialized kernel per slice (sizes may differ per slice) |
| **A_all_init_matrix** | cell {K×1} | 3D arrays [H×W×S] per kernel, unified via A1_matrix_unify_size (if use_matrix_format) |
| **init_kernel_centers** | array [K×2] | Centers [y,x] used (from IS01A for proliferation; from ref for block_manual) |

---

## All-Slice Decomposition Parameters (Block)

### Block: DA01A (Decompose-All-01-A)
### Function: `decomposeAllSlices.m` (wrapper); core: `MTSBD_synthetic_all_slice.m`

Standardized wrapper: uses organizeData/organizeParams for inputs/outputs. Algorithm parameters (initial_iteration, maxIT, lambda1, phase2_enable, lambda2, nrefine, kplus_factor, signflip_threshold, xpos, getbias, Xsolve_method) are inherited from params.slice (DS01A). Block-specific presets below.

| Parameter | Type | Units | Range | Default | Description |
|-----------|------|-------|-------|---------|-------------|
| **use_reference_init** | logical | - | {true, false} | true | Initialize X from reference slice (DS01A) results |
| **show_allslice_progress** | logical | - | {true, false} | true | Display optimization progress during all-slice decomposition |
| **maxIT_allslice** | integer | iterations | > 0 | 15 | Maximum outer iterations for all-slice decomposition |

#### Output (in `data.block`)

| Field | Type | Description |
|-------|------|-------------|
| **Aout_all** | cell {K×1} | Deconvolved kernels [H×W×S] per kernel |
| **Xout_all** | array [H×W×K] | Activation maps |
| **bout_all** | array [S×K] | Bias terms |
| **extras_all** | struct | Optimization details (phase1, phase2, runtime, etc.) |

**Requires:** IB01A (data.block.A_all_init_matrix), DS01A (data.slice for ref init and algo params), data.Y, data.X0.

---

## Visualization Parameters

### Block: VR01A (Visualize-Results-01-A)
### Function: `visualizeResults.m`

| Parameter | Type | Units | Range | Default | Description |
|-----------|------|-------|-------|---------|-------------|
| **visualize_all_slices** | logical | - | {true, false} | true | Visualize all slices (true) or reference only (false) |
| **save_figures** | logical | - | {true, false} | false | Save figure files to disk |

---

## Workspace Management Parameters

### Block: WS01A (Write-Save-01-A)

| Parameter | Type | Units | Range | Default | Description |
|-----------|------|-------|-------|---------|-------------|
| **save_workspace** | logical | - | {true, false} | true | Save complete MATLAB workspace |
| **save_results_only** | logical | - | {true, false} | false | Save only essential results (smaller file) |
| **results_path** | string | - | - | pwd | Directory path for saving results |

---

## Data Structures

### `data` Struct (Multi-dimensional Arrays)

Output from wrapper functions containing all data arrays:

```matlab
% After GD01A (Generate Data)
data.Y              % [H × W × S] - Normalized observation
data.A0             % {S × K} - Noisy kernels cell array (ground truth)
data.A0_noiseless   % {S × K} - Noiseless kernels cell array (ground truth)
data.X0             % [H × W × K] - Ground truth activations
data.Y_ref          % [H × W] - Reference slice observation
data.X0_ref         % [H × W × K] - Reference activations
data.A0_ref         % {1 × K} - Reference kernels (ground truth)

% After IN01A (Initialize Kernels) – in data.slice
data.slice.A_init   % {1 × K} - Initialized kernels for reference slice
data.slice.init_kernel_centers % [K × 2] - Centers (empty for manual)

% After DS01A (Decompose Reference Slice) – in data.slice
data.slice.A        % {1 × K} - Deconvolved kernels for reference slice
data.slice.X        % [H × W × K] - Activation maps for reference slice
data.slice.b        % [K × 1] - Bias terms for each kernel

% After IS01A (Find Isolated Points) – in data.slice
data.slice.kernel_centers       % [K × 2] - Most isolated points [y, x]
data.slice.A0_used              % {S × K} - Ground truth kernels (padded if needed)
data.slice.most_isolated_points % {1 × K} - Selected centers (cell)

% After IB01A (Initialize Block)
data.block.A_all_init           % {S × K} - Initialized kernels per slice
data.block.A_all_init_matrix   % {K × 1} - 3D arrays [H × W × S] per kernel
data.block.init_kernel_centers  % [K × 2] - Centers used for block init

% After DA01A (Decompose All Slices) – in data.block
data.block.Aout_all   % {K×1} - Deconvolved kernels [H×W×S]
data.block.Xout_all   % [H×W×K] - Activations
data.block.bout_all   % [S×K] - Bias terms
data.block.extras_all % Optimization details
```

**Note:** Block init (IB01A) and DA01A write to params.block and data.block; blocks before IB01A do not write to block.

Where:
- H = height in pixels
- W = width in pixels
- S = number of slices
- K = number of kernels

### `params` Struct (Parameters and Metadata)

Stored hierarchically as `params.synGen`, `params.slice`, `params.block`; `organizeParams(..., 'extract')` flattens to a single struct for use in functions.

**Hierarchical layout:**
- **params.synGen** – Data generation (SNR, N_obs, num_slices, kernel_sizes, LDoS_path, …).
- **params.slice** – Single-slice: init (init_method, init_window, kernel_selection_type, kernel_sizes_ref) and decomposition (initial_iteration, maxIT, lambda1, phase2_enable, …).
- **params.block** – Block init (IB01A): block_init_method, A1_matrix_unify_size, window_type_proliferation, etc. DA01A: maxIT_allslice, use_reference_init, show_allslice_progress. Written by IB01A and DA01A.

**Flat (after extract):**
```matlab
params.SNR, params.N_obs, params.num_slices, params.LDoS_path, params.num_kernels, params.ref_slice
params.init_method, params.init_window, params.initial_iteration, params.maxIT, params.lambda1
params.phase2_enable, params.xpos, params.Xsolve_method, ...
```

### `extras` Struct (Detailed Results)

Contains iteration history and detailed optimization information:

```matlab
extras.phase1.activation_metrics    % [IT × K] similarity over iterations
extras.phase1.kernel_quality_factors % [IT × K] quality over iterations
extras.phase1.Aout                  % Intermediate kernel results
extras.phase1.Xout                  % Intermediate activation results
extras.phase1.biter                 % Bias values per kernel

extras.phase2.*                     % (if phase2_enable = true)
extras.runtime                      % Total execution time
extras.normA                        % Kernel normalizations
```

---

## Parameter Value Guidelines

### Typical Working Ranges

| Parameter | Conservative | Moderate | Aggressive | Notes |
|-----------|--------------|----------|------------|-------|
| **SNR** | 5-10 | 3-5 | 1-3 | Lower = more challenging |
| **defect_density** | 0.001-0.005 | 0.01-0.03 | 0.05-0.1 | Higher = more crowded |
| **lambda1** | 1e-1 to 5e-2 | 3e-2 to 1e-2 | 5e-3 to 1e-3 | Lower = less regularization |
| **maxIT** | 10-15 | 20-30 | 50+ | More iterations = slower but potentially better |
| **observation_resolution** | 3 | 4-5 | 6+ | Higher = more detail but slower |

### Parameter Dependencies

- **lambda1 ∝ 1/SNR**: Lower SNR requires higher regularization
- **maxIT ∝ num_kernels**: More kernels may need more iterations
- **kplus_factor × kernel_size**: Should be < min(H,W)/4 to avoid boundary issues
- **defect_density × N_obs²**: Total defects ≈ this product

---

## Workspace Management Parameters

### Block: WS01A (Workspace-Save-01-A)
### Functions: `saveDataset.m`, `saveRun.m`, `autoSave.m`

**Design**: Always executes when reached - no enable flag needed. Comment out block to skip.

**Note**: `saveWorkspace.m` has been removed. Use `saveDataset.m` (pre-run) or `saveRun.m` (post-run) instead. `autoSave.m` automatically selects the appropriate function based on workflow phase.

**Automatic Saving**: The workflow uses structured saves:
- **Pre-run phase**: `saveDataset.m` saves to `project/<dataset>/<dataset>.mat`
- **Post-run phase**: `saveRun.m` saves to `project/<dataset>/runXX/runXX.mat`
- **Auto-save**: `autoSave.m` automatically calls the correct function

See `saveDataset.m` and `saveRun.m` documentation for details.

### Block: WL01A (Workspace-Load-01-A)
### Function: `loadWorkspace.m`

**Design**: Always executes when reached - no enable flag needed. Comment out block to skip.

**Memory Protection**: Always prompts user before overwriting current workspace in memory.

| Parameter | Type | Units | Range | Default | Description |
|-----------|------|-------|-------|---------|-------------|
| **workspace_dir** | string | - | valid path | '' (UI selection) | Directory containing workspace. Empty = interactive folder picker |
| **verify_log_file** | logical | - | {true, false} | true | Verify log file naming convention matches |

#### Loading Behavior

1. **Directory Selection**: UI dialog if `workspace_dir` empty
2. **File Detection**: Automatically finds .mat file in directory
3. **Naming Validation**: Verifies `<name>.mat` and `<name>_LOGfile.txt` match
4. **Collision Handling**: If current workspace has same name, prompts user:
   - Save current workspace first (recommended)
   - Rename and save current workspace
   - Delete current workspace (with confirmation)
   - Cancel load operation

#### Naming Convention (Enforced)

**Required Structure**:
```
workspace_directory/
  ├─ experiment_001.mat          ← Workspace file
  └─ experiment_001_LOGfile.txt  ← Log file (MUST match name)
```

**Error**: If log file name doesn't match, load will fail with clear error message

---

## Version History

- **v1.1** (2025-10-27): Updated IN01A block
  - Added `show_ground_truth` and `show_initialized` parameters
  - Updated to use `initializeKernelsRef` wrapper function
  - Clarified kernel_sizes_ref output in params

- **v1.0** (2025-10-27): Initial glossary
  - Documented GD01A (Data Generation) parameters
  - Documented IN01A (Kernel Initialization) parameters
  - Documented DS01A (Decomposition) parameters
  - Documented IS01A, IB01A (block init), DA01A, VR01A, WS01A parameters
  - Defined data structures (data, params, extras)

---

## Notes

- All parameters are case-sensitive
- Parameters marked as *required* must be provided (or user will be prompted)
- Default values are suggestions; optimal values depend on specific dataset
- Units marked as "-" are dimensionless
- Ranges shown as (a,b) exclude endpoints; [a,b] include endpoints

---

**Maintained by**: Dong Chen  
**For questions**: See function documentation or README.md

