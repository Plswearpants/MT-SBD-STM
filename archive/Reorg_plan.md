
  # MT-SBD-STM Codebase Reorganization Plan

  ## Goal
  Reorganize ~200 MATLAB files from a flat/disorganized layout into a clean top-down hierarchy: **user-facing scripts → algorithm API → workflow wrappers → core solvers → internal helpers**. A new user
  should understand the codebase by looking at the top-level directory.

  ## Proposed Structure

  ```
  MT-SBD-STM/
  ├── init_sbd.m                 (ONLY user-facing .m at root)
  ├── startup.m                  (MATLAB auto-startup)
  ├── README.md / CLAUDE.md / .gitignore
  │
  ├── scripts/                   [LAYER 1: User workflows - run these]
  │   ├── run_synthetic_data.m
  │   ├── MTSBD_block_synthetic.m
  │   ├── MTSBD_block_realdata.m
  │   └── postprocessing.m  (+ 2 illustration scripts)
  │
  ├── examples/                  [LAYER 1: Demos and tutorials]
  │   ├── simple_SBD_example.m
  │   ├── (existing examples stay in place)
  │   ├── phase_space/
  │   ├── data/                  (material-specific: Ag111, ZrSiTe, etc.)
  │   └── example_data/          (tracked .mat reference data)
  │
  ├── src/                       [LAYERS 2-5: All internal source code]
  │   ├── algorithms/            [LAYER 2: Public API functions]
  │   │   ├── MT_SBD.m
  │   │   ├── SBD.m
  │   │   ├── MTSBD_synthetic.m
  │   │   └── MTSBD_*.m         (all 7 algorithm variants)
  │   │
  │   ├── workflow/              [LAYER 3: High-level wrappers]
  │   │   ├── generateSyntheticData.m
  │   │   ├── decomposeReferenceSlice.m
  │   │   ├── findIsolatedPoints.m
  │   │   ├── saveDataset.m / saveRun.m / loadWorkspace.m
  │   │   └── (13 files from Dong_func/wrapper/)
  │   │
  │   ├── solvers/               [LAYER 4: Optimization solvers]
  │   │   ├── Asolve_Manopt*.m   (5 variants)
  │   │   ├── Xsolve_FISTA*.m    (5 variants)
  │   │   ├── Xsolve_pdNCG.m
  │   │   └── trustregions.m
  │   │
  │   ├── helpers/               [LAYER 5: Low-level math primitives]
  │   │   ├── convfft2.m / convfft3.m / cconvfft2.m
  │   │   ├── H_function.m / Hxx_function.m
  │   │   ├── huber.m / shrink.m / extend.m
  │   │   └── proj2oblique.m / seq_crosscorr_regularizer.m
  │   │
  │   ├── initialization/        (kernel init functions)
  │   ├── metrics/               (quality metrics & evaluation)
  │   ├── generation/            (synthetic data generation)
  │   ├── preprocessing/         (streak removal, normalization, masking)
  │   ├── io/                    (STM file loading, parallel result loading)
  │   ├── visualization/         (all plotting and display functions)
  │   └── utils/                 (general-purpose utilities, logging)
  │
  ├── config/                    [Configuration - centralized]
  │   ├── default_config_settings.m
  │   └── *.mat config files     (generated at runtime)
  │
  ├── tests/                     [All test/experimental scripts]
  │   ├── SBD_test*.m            (5 test scripts from root)
  │   ├── SBD_with_proximity_script.m
  │   ├── Asolve_Manopt_tunable_test.m
  │   └── Xsolve_FISTA_test.m
  │
  ├── lib/                       [Third-party code]
  │   ├── imshow3D/
  │   ├── colormap/
  │   └── mat2im/
  │
  ├── docs/                      (documentation)
  ├── template/                  (script templates)
  ├── archive/                   (historical/deprecated: multi_kernel_test.m, etc.)
  ├── tools/                     (developer utilities: function_path_analyzer.m)
  ├── projects/                  (runtime output - gitignored)
  └── session/                   (runtime output - gitignored)
  ```

  ---

  ## Implementation Steps

  Each step produces a single git commit. Steps 1-7 move files; Step 8 updates path infrastructure; Steps 9-10 fix references and docs.

  ### Step 0: Create safety branch ✅ DONE                                                                                                                                                                         
  - Already on `ClaudeCode` branch — serves as our working branch, with `main` as rollback.

  ### Step 1: Create directory skeleton
  - Create all empty directories: `src/algorithms/`, `src/solvers/`, `src/helpers/`, `src/workflow/`, `src/initialization/`, `src/metrics/`, `src/generation/`, `src/preprocessing/`, `src/io/`,
  `src/visualization/`, `src/utils/`, `tests/`, `lib/imshow3D/`, `lib/colormap/`, `lib/mat2im/`, `archive/`, `tools/`
  - No files move yet.
  - **Commit**: `chore: create new directory skeleton`

  ### Step 2: Move third-party libraries to `lib/`
  - `Dong_func/imshow3D.m/*` → `lib/imshow3D/`
  - `colormap/*` → `lib/colormap/`
  - `utils/mat2im/*` → `lib/mat2im/`
  - **Commit**: `refactor: move third-party libraries to lib/`

  ### Step 3: Move core helpers to `src/helpers/`
  - All 10 files from `core/helpers/` → `src/helpers/`
  - **Critical path fix**: `H_function.m` line 3 loads `'/../../config/Hfunction_config.mat'` — from `src/helpers/` this resolves to root `config/`, same as before. No change needed.
  - **Commit**: `refactor: move core helpers to src/helpers/`

  ### Step 4: Move solvers to `src/solvers/`
  - 11 solver files from `core/` → `src/solvers/`
  - 2 test files from `core/` → `tests/` (`Asolve_Manopt_tunable_test.m`, `Xsolve_FISTA_test.m`)
  - Move stray tunable config .mat files from root and `examples/` → `config/`

  **Critical path fixes in 14 files** (the main risk of this reorganization):

  | File | Old pattern | New pattern |
  |------|------------|-------------|
  | `Asolve_Manopt.m` | `'/../config/Asolve_config.mat'` | `'/../../config/Asolve_config.mat'` |
  | `Asolve_Manopt_ALL.m` | `'/../config/Asolve_config.mat'` | `'/../../config/Asolve_config.mat'` |
  | `Asolve_Manopt_multikernel.m` | `'/../config/Asolve_config.mat'` | `'/../../config/Asolve_config.mat'` |
  | `Asolve_Manopt_tunable.m` | `'/../examples/Asolve_config_tunable.mat'` | `'/../../config/Asolve_config_tunable.mat'` |
  | `Xsolve_FISTA.m` | `[fpath '/../config/...']` | `[fpath '/../../config/...']` |
  | `Xsolve_FISTA_tunable.m` | `'/../examples/Xsolve_config_tunable.mat'` | `'/../../config/Xsolve_config_tunable.mat'` |
  | `Xsolve_FISTA_ALL.m` | `[fpath '/../config/...']` | `[fpath '/../../config/...']` |
  | `Xsolve_FISTA_multikernel.m` | `[fpath '/../config/...']` | `[fpath '/../../config/...']` |
  | `Xsolve_FISTA_parallel.m` | `[fpath '/../config/...']` | `[fpath '/../../config/...']` |
  | `Xsolve_pdNCG.m` | `[fpath '/../config/...']` | `[fpath '/../../config/...']` |

  Also in all Xsolve files: remove `addpath([fpath '/helpers'])` lines (unnecessary once `init_sbd.m` adds `src/` recursively).

  - **Commit**: `refactor: move solvers to src/solvers/ and fix config paths`

  ### Step 5: Move algorithm functions to `src/algorithms/`
  - 7 files from root → `src/algorithms/`: `MT_SBD.m`, `SBD.m`, `MTSBD_synthetic.m`, `MTSBD_synthetic_Xregulated.m`, `MTSBD_synthetic_all_slice.m`, `MTSBD_synthetic_Xregulated_all_slices.m`,
  `MTSBD_all_slice.m`
  - No path fixes needed (these use function names, not relative paths).
  - **Commit**: `refactor: move algorithm functions to src/algorithms/`

  ### Step 6: Decompose `Dong_func/` into `src/` modules
  This is the largest step — moves ~95 files from the `Dong_func/` catch-all into proper modules.

  **Moves by target:**
  - `Dong_func/wrapper/` (13 files) → `src/workflow/`
  - `Dong_func/basic/` (5 files) → `src/utils/`
  - `Dong_func/data_preprocessing/` (2 files) → `src/preprocessing/`
  - Initialization functions (5) → `src/initialization/`
  - Quality metrics (13) → `src/metrics/`
  - Data generation (6) → `src/generation/`
  - Preprocessing/streaks (11) → `src/preprocessing/`
  - STM I/O functions (6) → `src/io/`
  - Visualization functions (13) → `src/visualization/`
  - Remaining utilities (~20) → `src/utils/`
  - Delete empty `Dong_func/` tree

  Note: `Dong_func/` was never added to the MATLAB path by `init_sbd.m`. Moving its contents into `src/` (which will be on the path) actually **fixes** this problem.

  - **Commit**: `refactor: decompose Dong_func/ into src/ modules`

  ### Step 7: Move remaining root files and old `utils/`
  - 5 test scripts from root → `tests/` (`SBD_test*.m`)
  - `SBD_with_proximity_script.m` → `tests/`
  - `plot_activations.m` → `src/visualization/`
  - `streak_correction.m` → `src/preprocessing/`
  - `function_path_analyzer.m` → `tools/`
  - `plotting/streaknoise.m` → `src/visualization/`
  - Old `utils/` (16 files: showims.m, Get_Colors.m, etc.) → `src/utils/`
  - `historical/*` → `archive/`
  - Move stray root `.mat` data files to `examples/data/` or gitignore
  - Delete empty old directories: `core/`, `utils/`, `colormap/`, `plotting/`, `historical/`, `Dong_func/`

  After this step, root contains ONLY: `init_sbd.m`, `startup.m`, `README.md`, `CLAUDE.md`, `.gitignore`, and directories.

  - **Commit**: `refactor: move remaining root files and clean up`

  ### Step 8: Update `init_sbd.m` and `startup.m`
  Rewrite `init_sbd.m` to add the new directory tree:
  ```matlab
  % Old: addpath(genpath([fp 'core'])), addpath(genpath([fp 'utils'])), addpath(genpath([fp 'config']))
  % New:
  addpath(genpath([fp 'src']));      % All internal source code
  addpath(genpath([fp 'config']));   % Configuration
  addpath(genpath([fp 'lib']));      % Third-party libraries
  ```

  Update `startup.m`: replace `addpath(genpath(fullfile(pwd, 'colormap')))` with `addpath(genpath(fullfile(pwd, 'lib')))` (or make it a no-op since `init_sbd` handles everything).

  - **Commit**: `refactor: update init_sbd.m for new directory structure`

  ### Step 9: Fix `run('../init_sbd')` references
  16 files reference `run('../init_sbd')`. Since `scripts/` and `examples/` remain at depth 1 from root, most references are already correct. Fix:
  - `examples/phase_space/phasespace_dataset.m`: `run('../init_sbd')` → `run('../../init_sbd')` (existing bug — it's at depth 2)

  - **Commit**: `fix: correct init_sbd path references in scripts`

  ### Step 10: Update documentation
  - Rewrite `CLAUDE.md` to reflect new structure
  - Update `README.md` with new directory layout
  - Update `docs/FOLDER_STRUCTURE_ANALYSIS.md`
  - Update `.gitignore` if needed for new paths

  - **Commit**: `docs: update documentation for new directory structure`

  ---

  ## Critical Files to Modify

  | File | Change |
  |------|--------|
  | `init_sbd.m` | Rewrite path setup for `src/`, `config/`, `lib/` |
  | `startup.m` | Update colormap path to `lib/colormap/` |
  | `core/Asolve_Manopt.m` | `/../config/` → `/../../config/` |
  | `core/Asolve_Manopt_ALL.m` | Same |
  | `core/Asolve_Manopt_multikernel.m` | Same |
  | `core/Asolve_Manopt_tunable.m` | `/../examples/` → `/../../config/` |
  | `core/Xsolve_FISTA.m` | `/../config/` → `/../../config/`, remove `addpath helpers` |
  | `core/Xsolve_FISTA_tunable.m` | `/../examples/` → `/../../config/`, remove `addpath helpers` |
  | `core/Xsolve_FISTA_ALL.m` | Same pattern |
  | `core/Xsolve_FISTA_multikernel.m` | Same pattern |
  | `core/Xsolve_FISTA_parallel.m` | Same pattern |
  | `core/Xsolve_pdNCG.m` | Same pattern |
  | `examples/phase_space/phasespace_dataset.m` | Fix `run('../init_sbd')` → `run('../../init_sbd')` |
  | `CLAUDE.md` | Rewrite for new structure |

  ---

  ## Verification Plan

  After each step, verify the reorganization hasn't broken anything:

  1. **After Step 8** (path updates): Run `init_sbd('verbose')` from root — should print success and have no errors
  2. **Spot-check function resolution**:
  - `which MT_SBD` → `src/algorithms/MT_SBD.m`
  - `which Asolve_Manopt_tunable` → `src/solvers/Asolve_Manopt_tunable.m`
  - `which convfft2` → `src/helpers/convfft2.m`
  - `which generateSyntheticData` → `src/workflow/generateSyntheticData.m`
  - `which computeQualityMetrics` → `src/metrics/computeQualityMetrics.m`
  - `which showims` → `src/utils/showims.m`
  - `which invgray` → `lib/colormap/invgray.m`
  3. **End-to-end**: Run `scripts/run_synthetic_data.m` through the first few blocks
  4. **Example**: Run `examples/simple_SBD_example.m` end-to-end
  5. **Root cleanliness**: `ls *.m` at root returns only `init_sbd.m` and `startup.m`

  ---

  ## Risk Mitigations

  - **Git branch**: All work on `reorganize/codebase-cleanup`, main branch untouched
  - **Stepwise commits**: Each step is independently committable and reversible
  - **`fileparts` paths**: All 14 solver files with self-relative config paths are explicitly updated in Step 4
  - **`Dong_func/` not on path**: Moving its contents to `src/` (which IS on the path) fixes an existing problem
  - **Config .mat duplication**: Centralizing all config files in `config/` eliminates the current scatter across root, `core/`, `examples/`, and `scripts/`
  - **Function name conflicts**: No duplicates exist across the moved files (verified)