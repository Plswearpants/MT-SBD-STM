# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

MT-SBD-STM implements **Multi-Type Sparse Blind Deconvolution** for Scanning Tunneling Microscopy (STM) images in MATLAB. It extends single-kernel SBD (Cheung et al., 2020) to simultaneously deconvolve multiple defect types from STM observations using the Riemannian Trust-Region Method (RTRM).

The core problem: given an observation Y containing multiple overlapping defect signatures, decompose it into kernels A (defect patterns) and activations X (spatial positions) such that Y ≈ Σ(A_n * X_n).

## Setup and Running

**Prerequisites**: MATLAB with Signal Processing and Image Processing Toolboxes. The [Manopt](http://www.manopt.org) library must be installed and on the MATLAB path.

**Initialization**: Every session must call `init_sbd()` from the repo root. This verifies Manopt is available, adds `src/`, `config/`, `lib/` to the path, and loads default solver configuration.

**Primary workflows**:
- `scripts/run_synthetic_data.m` — Standardized synthetic data pipeline (block-based, with logging)
- `scripts/MTSBD_block_synthetic.m` — Block-based synthetic processing (older)
- `scripts/MTSBD_block_realdata.m` — Real experimental STM data processing
- `examples/MT_SBD_STM_experimental_data.m` — Real data from .3ds files

**Running a workflow**: Open MATLAB in the repo root, ensure Manopt is on the path, then run `init_sbd` followed by the desired script. Scripts use block-based structure — run individual blocks (MATLAB sections) with Ctrl+Enter.

**Testing**: No formal test framework. Validation is done via test scripts in `tests/`:
- `tests/SBD_test*.m` — Comprehensive algorithm tests
- `examples/MT_SBD_STM_synthetic_testing.m` — Synthetic data validation

## Directory Structure

```
MT-SBD-STM/
├── init_sbd.m                 (Entry point - adds paths)
├── startup.m                  (MATLAB auto-startup)
│
├── scripts/                   [User workflows - run these]
│   ├── run_synthetic_data.m
│   ├── MTSBD_block_synthetic.m
│   └── MTSBD_block_realdata.m
│
├── examples/                  [Demos, tutorials, experimental data]
│   ├── simple_SBD_example.m
│   ├── phase_space/           (parameter space analysis)
│   └── data/                  (material-specific: Ag111, ZrSiTe, etc.)
│
├── src/                       [All internal source code]
│   ├── algorithms/            [Public API: MT_SBD.m, SBD.m, MTSBD_*.m]
│   ├── solvers/               [Optimization: Asolve_*, Xsolve_*, trustregions]
│   ├── helpers/               [Math primitives: convfft2, H_function, shrink]
│   ├── workflow/              [High-level wrappers: generateSyntheticData, save*]
│   ├── initialization/        [Kernel initialization functions]
│   ├── metrics/               [Quality metrics & evaluation]
│   ├── generation/            [Synthetic data generation]
│   ├── preprocessing/         [Streak removal, normalization, masking]
│   ├── io/                    [STM file loading, parallel result loading]
│   ├── visualization/         [Plotting and display functions]
│   └── utils/                 [General utilities, logging]
│
├── config/                    [Configuration files]
├── lib/                       [Third-party: imshow3D, colormap, mat2im]
├── tests/                     [Test scripts]
├── archive/                   [Historical/deprecated code]
└── tools/                     [Developer utilities]
```

## Architecture

### Core Algorithm (`src/algorithms/MT_SBD.m`, `src/algorithms/SBD.m`)

Two-phase alternating optimization:

**Phase I** — Alternates between kernel solver (Asolve) and activation solver (Xsolve) for `maxIT` iterations. Uses a demixing strategy: for each kernel n, compute the residual Y - Σ(other kernels) and solve for A_n and X_n on that residual. Kernels are processed in order of descending variance.

**Phase II** (optional) — Sphere lifting: pads kernels by `kplus` pixels to allow spatial refinement, then performs lambda continuation from `lambda1` → `lambda2` over `nrefine` steps. Includes unshifting via circshift to find optimal spatial alignment.

### Solver Components (`src/solvers/`)

- **`Asolve_Manopt*.m`** — Kernel sub-problem solver using Manopt's trust-region optimizer on the unit sphere manifold. Tunable variant (`_tunable`) accepts solver config at runtime.
- **`Xsolve_FISTA*.m`** — Activation sub-problem using Fast Iterative Shrinkage-Thresholding Algorithm. Supports positivity constraints and bias extraction.
- **`Xsolve_pdNCG.m`** — Alternative primal-dual nonlinear conjugate gradient solver for activations.
- **`trustregions.m`** — Core Riemannian trust-region solver (from Manopt, customized).

### Helper Functions (`src/helpers/`)

- FFT-based convolution: `convfft2.m`, `convfft3.m`, `cconvfft2.m`
- Hessian computation: `H_function.m`, `Hxx_function.m`
- Proximal operators: `shrink.m`, `huber.m`

### Source Modules (`src/`)

- **`workflow/`** — High-level wrappers (e.g., `createProjectStructure.m`, `saveDataset.m`, `saveRun.m`, `autoSave.m`, `loadWorkspace.m`, `generateSyntheticData.m`)
- **`preprocessing/`** — Streak detection/correction, normalization, masking
- **`metrics/`** — Quality metrics computation and evaluation
- **`generation/`** — Synthetic data generation utilities
- **`visualization/`** — Plotting and display functions
- **`utils/`** — General utilities including `logUsedBlocks.m` (logging system)

### Configuration (`config/`)

- `default_config_settings.m` — Sets Manopt (verbosity, tolgradnorm, maxiter) and FISTA (MAXIT, EPSILON) defaults
- `.mat` config files — Binary solver configuration (Asolve, Xsolve, Hfunction)
- Runtime config updates via `update_config.m`

### Key Data Structures

```
data.Y              — [H × W × S] observation (S = energy slices)
data.A0             — {S × K} ground truth kernels (cell array)
data.X0             — [H × W × K] ground truth activations
params              — Struct: SNR, N_obs, lambda1/2, phase2, xpos, getbias, Xsolve, etc.
extras              — Struct: per-phase quality_metrics, residuals, Aout, Xout, runtime
```

## Script Block Convention

Scripts in `scripts/` use a block-based architecture with `ABXXZ` identifiers:
- **A** = Category (G=Generate, I=Initialize, D=Decompose, V=Visualize, W=Write)
- **B** = Subcategory (D=Data, N=kNernel, S=Slice, P=Proliferation, R=Results)
- **XX** = Sequential number, **Z** = Variant letter

Each block has a PRESETS section (user-configurable), a DO NOT EDIT section (logic), and ends with `clearvars` to prevent namespace pollution. Blocks are logged via `logUsedBlocks()`.

## Project Output Structure

Results are stored in `projects/synthetic_<timestamp>/`:
```
projects/synthetic_<timestamp>/
  ├─ auto.mat + auto_LOGfile.txt           (dataset + log)
  ├─ run_auto_01/run_auto_01.mat + log     (algorithm run results)
  └─ run_auto_02/...
```

## Key Conventions

- All convolutions use FFT-based operations (`convfft2.m`, `convfft3.m`), not spatial convolution
- Kernels are constrained to the unit sphere (Riemannian manifold optimization)
- Lambda parameters can be scalar (all kernels) or vector (per-kernel)
- The `_tunable` suffix on solver functions indicates runtime-configurable variants
- `.mat` files and `projects/`, `session/` directories are gitignored; example LDoS data files are explicitly tracked
- Synthetic workflows include ground truth comparison with quality metrics (activation similarity, kernel quality factors)
