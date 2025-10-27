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
5. [Proliferation Parameters](#proliferation-parameters)
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
### Function: `initialize_kernels.m`

| Parameter | Type | Units | Range | Default | Description |
|-----------|------|-------|-------|---------|-------------|
| **kernel_selection_type** | string | - | {'selected', 'random'} | 'selected' | Method for initial kernel selection. 'selected' = interactive, 'random' = random initialization |
| **window_type** | cell/string | - | see below | {'gaussian', 2.5} | Window function applied to kernels to reduce edge effects |

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
### Functions: `MTSBD_synthetic.m`, `Asolve_Manopt_tunable.m`, `Xsolve_FISTA_tunable.m`

#### Algorithm Selection

| Parameter | Type | Units | Range | Default | Description |
|-----------|------|-------|-------|---------|-------------|
| **use_Xregulated** | logical | - | {true, false} | false | Use X-regularized version with cross-correlation penalty |

#### Phase I Settings (Initial Decomposition)

| Parameter | Type | Units | Range | Default | Description |
|-----------|------|-------|-------|---------|-------------|
| **initial_iteration** | integer | iterations | ≥ 1 | 1 | Number of inner solver iterations at start. Increases gradually |
| **maxIT** | integer | iterations | > 0 | 15 | Maximum number of outer alternating iterations |
| **lambda1** | vector [K×1] | - | > 0 | [3e-2, ...] | L1 regularization parameter for Phase I (per kernel). Controls sparsity |

#### Phase II Settings (Refinement - Optional)

| Parameter | Type | Units | Range | Default | Description |
|-----------|------|-------|-------|---------|-------------|
| **phase2_enable** | logical | - | {true, false} | false | Enable Phase II refinement with sphere lifting |
| **lambda2** | vector [K×1] | - | > 0 | [1e-2, ...] | Final L1 regularization for Phase II (per kernel) |
| **nrefine** | integer | steps | ≥ 1 | 5 | Number of refinement steps. Lambda decreases from lambda1 to lambda2 |
| **kplus_factor** | float | - | > 0 | 0.5 | Sphere lifting padding factor: `kplus = kplus_factor × kernel_size` |

#### Algorithm Parameters

| Parameter | Type | Units | Range | Default | Description |
|-----------|------|-------|-------|---------|-------------|
| **signflip_threshold** | float | - | [0, 1] | 0.2 | Threshold for detecting sign flips in optimization |
| **xpos** | logical | - | {true, false} | true | Enforce non-negative constraint on activations X |
| **getbias** | logical | - | {true, false} | true | Extract constant bias term from observation |
| **Xsolve_method** | string | - | {'FISTA', 'pdNCG'} | 'FISTA' | Solver for activation map X. FISTA = Fast Iterative Shrinkage-Thresholding |
| **gamma_crosscorr** | float | - | > 0 | 5e-2 | Cross-correlation regularization parameter (only if use_Xregulated = true) |
| **use_xinit** | struct/[] | - | - | [] | Initial guess for X. Empty = start from scratch |

#### Quality Metrics (Output)

| Parameter | Type | Units | Description |
|-----------|------|-------|-------------|
| **activation_metrics** | matrix [IT×K] | - | Activation similarity at each iteration for each kernel |
| **kernel_quality_factors** | matrix [IT×K] | - | Kernel quality score at each iteration for each kernel |

---

## Isolation Analysis Parameters

### Block: IS01A (Isolation-Selection-01-A)
### Function: TBD (to be wrapped)

| Parameter | Type | Units | Range | Default | Description |
|-----------|------|-------|-------|---------|-------------|
| **isolation_threshold_factor** | float | - | > 1 | 10 | Defect detection threshold: `threshold = max_activation / factor` |
| **target_kernel_size_type** | string | - | see below | 'kernel_sizes_all' | Method for determining target kernel sizes across slices |

#### Target Kernel Size Options:
- `'ref_kernel_sizes'` - Use reference slice kernel sizes for all slices
- `'kernel_sizes_cap'` - Use maximum kernel size across all slices
- `'kernel_sizes_all'` - Use slice-specific kernel sizes

#### Derived/Output Parameters

| Parameter | Type | Units | Description |
|-----------|------|-------|-------------|
| **most_isolated_points** | cell {1×K} | pixels | Coordinates [y, x] of most isolated defect for each kernel |
| **isolation_scores** | cell {1×K} | pixels² | Squared distance to nearest other-kernel defect |
| **defect_positions** | cell {1×K} | pixels | All detected defect positions [N×2] for each kernel |
| **num_defects** | vector [1×K] | - | Number of defects detected per kernel |

---

## Proliferation Parameters

### Block: IP01A (Initialize-Proliferation-01-A)
### Function: `initialize_kernels_proliferation.m`

| Parameter | Type | Units | Range | Default | Description |
|-----------|------|-------|-------|---------|-------------|
| **window_type_proliferation** | cell | - | see IN01A | {'gaussian', 2.5} | Window function for proliferated kernels |
| **interactive_size_adjust** | logical | - | {true, false} | false | Allow interactive kernel size adjustment during proliferation |
| **use_matrix_format** | logical | - | {true, false} | true | Output format: true = 3D matrix [H×W×S], false = cell array {S×K} |

---

## All-Slice Decomposition Parameters

### Block: DA01A (Decompose-All-01-A)
### Function: `MTSBD_synthetic_all_slice.m`

| Parameter | Type | Units | Range | Default | Description |
|-----------|------|-------|-------|---------|-------------|
| **use_Xregulated_allslice** | logical | - | {true, false} | false | Use X-regularized version for all-slice processing |
| **use_reference_init** | logical | - | {true, false} | true | Initialize X from reference slice results |
| **show_allslice_progress** | logical | - | {true, false} | true | Display optimization progress during all-slice decomposition |
| **maxIT_allslice** | integer | iterations | > 0 | 15 | Maximum iterations for all-slice decomposition |

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
data.Y              % [H × W × S] - Normalized observation
data.A0             % {S × K} - Noisy kernels cell array
data.A0_noiseless   % {S × K} - Noiseless kernels cell array
data.X0             % [H × W × K] - Ground truth activations
data.Y_ref          % [H × W] - Reference slice observation
data.X0_ref         % [H × W × K] - Reference activations
data.A0_ref         % {1 × K} - Reference kernels
```

Where:
- H = height in pixels
- W = width in pixels
- S = number of slices
- K = number of kernels

### `params` Struct (Parameters and Metadata)

Contains all scalar parameters, strings, and small arrays:

```matlab
% Input parameters
params.SNR                      % Signal-to-noise ratio
params.N_obs                    % Observation lattice size
params.observation_resolution   % Pixels per lattice site
params.defect_density           % Surface defect density
params.num_slices               % Number of energy slices
params.LDoS_path                % Path to LDoS data

% Derived parameters
params.num_kernels              % Number of kernels
params.ref_slice                % Reference slice index
params.kernel_sizes             % [S × K × 2] kernel dimensions

% Algorithm parameters
params.lambda1                  % [K × 1] regularization
params.phase2                   % Boolean for Phase II
params.xpos                     % Boolean for positivity
params.Xsolve                   % Solver type string
% ... etc
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

## Version History

- **v1.0** (2025-10-27): Initial glossary
  - Documented GD01A (Data Generation) parameters
  - Documented IN01A (Kernel Initialization) parameters
  - Documented DS01A (Decomposition) parameters
  - Documented IS01A, IP01A, DA01A, VR01A, WS01A parameters
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

