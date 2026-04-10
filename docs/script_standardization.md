## MT-SBD-STM script standardization (draft)

This file defines a **minimal, high-reliability** standard for scripts in `scripts/` (e.g. `run_synthetic_data.m`, `run_real_data.m`). The goal is to make scripts **repeatable**, **auditable**, and easy to refactor without changing scientific intent.

### Goals (non-negotiable)

- **Reproducible**: all user choices are captured into `params`/`cfg` and saved, so a user can restart from any checkpoint with no missing inputs.
- **Traceable**: each executed block is logged with `logUsedBlocks(...)`.
- **Minimal surface area**: scripts orchestrate; heavy logic lives in functions (wrappers) under `Dong_func/` or `core/`/`utils/`.
- **Portable paths**: no hard-coded absolute paths unless explicitly required; prefer `fullfile(...)` + repo-relative discovery.

### Project-first UI model (checkpoint tree)

The **primary UX** is “choose a project, then choose a checkpoint node to load”.

- **Project**: a folder that contains a coherent analysis lifecycle and its saved artifacts.
- **Checkpoint**: a saved state that a user can load and continue from (node in the tree).
- **Branching**: continuing work from a checkpoint creates a *new child* checkpoint (or a new run folder) rather than overwriting the parent. This supports “what-if” exploration without losing history.

#### Proposed checkpoint nodes (simple + intuitive)

Use a small, stable set of checkpoint types (synthetic/real can share most of these):

- **raw**: raw loaded measurement (e.g. `.3ds`) + minimal metadata.
- **preprocess**: preprocessing outputs (cropped/masked/normalized data) + preprocessing params.
- **pre-run**: initialization artifacts before main solve (reference slice choice, initial kernels, initial activations, noise estimates).
- **run**: algorithm outputs (kernels/activations/bias/extras) + solver params.
- **post-run**: derived analyses (quality metrics, reconstructions, summary figures/tables).

The script standardization goal is to make each node **loadable** and **resumable**.

### Project & checkpoint storage model

- **Exactly one project per physical dataset** (per user workflow):
  - All branches/checkpoints that depend on a given raw `.3ds` live under the same project root.
- **Exactly two raw copies**:
  - Copy 1: the original lab/raw-data file (immutable, never altered).
  - Copy 2: a single **project-local copy** (e.g. under `projects/real_<timestamp>_<basename>/raw/`).
  - All checkpoints in the project reference this project-local raw; they do **not** re-save the raw.
- **Incremental checkpoints**:
  - Each checkpoint file stores only:
    - the **new stage’s results + params** (e.g. preprocessing outputs and parameters),
    - a pointer to the **project-local raw**,
    - a pointer to its **parent node**.
  - Minimal required metadata per node:
    - `meta.project_root` (project folder),
    - `meta.raw_path_project` (path to project-local raw),
    - `meta.node_id` (unique within project),
    - `meta.parent_id` (empty only for the root / raw node),
    - `meta.stage` (`raw | preprocess | pre-run | run | post-run`),
    - `meta.timestamp` (creation time),
    - `meta.label` (short human-readable title),
    - optional `meta.abstract` (1–3 line summary of what changed vs parent).

### Canonical script layout

1. **Section 0: Path initialization**
   - `clc; clear; close all;` (optional, but default for trunk scripts)
   - Run `init_sbd` from repo root:

```matlab
% scripts/<name>.m
clc; clear; close all;

repo_root = fileparts(fileparts(mfilename('fullpath'))); % .../<repo>/
init_path = fullfile(repo_root, 'init_sbd.m');
if exist(init_path, 'file')
    run(init_path);
else
    error('init_sbd.m not found. Run from repo or keep folder structure intact.');
end
```

2. **Block structure (repeatable pattern)**
   - Each block starts with a banner and a **block ID**.
   - Each block has a **PRESETS** region (only user-editable lines).
   - The rest is **DO NOT EDIT BELOW** (calls wrappers, logs, saves).

#### Block banner template

```matlab
%% =========================================================================
%% GD01A: Generate-Data-01-A; One-line summary
%  =========================================================================
%  Dependencies: <function1.m>, <function2.m>
%
% -------------------------------------------------------------------------
% PRESETS: User-configurable parameters
% -------------------------------------------------------------------------
% params.<group>.<name> = <value>;
%
%%%%%%%%%%%%%%%%%% DO NOT EDIT BELOW %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOGcomment = logUsedBlocks(log.path, log.file, "GD01A", LOGcomment, 0/1);
% [data, params] = wrapperFunction(log, data, params, ...);
```

### Data model (script-level variables)

- **`params`**: user-facing parameters, grouped by pipeline stage (`params.synGen`, `params.mcsbd_slice`, `params.mcsbd_block`, ...).
- **`data`**: observations and algorithm outputs (`data.synGen.Y`, `data.mcsbd_slice.*`, `data.mcsbd_block.*`, ...).
- **`meta`**: filesystem metadata (project/run folder names and paths; see `createProjectStructure`, `saveDataset`, `saveRun`).
- **`log`**: logging target for `logUsedBlocks` (`log.path`, `log.file`).

### Logging requirements (reliability)

- Every meaningful block must call `logUsedBlocks(...)` at **block start**.
- A session must create a **fresh** log file (avoid log contamination):
  - Use `initialize=1` only once per dataset/run session.
  - Then use `initialize=0` for all subsequent entries.
- Log comments should record **only** the critical presets and any interactive outcomes (e.g. selected file, chosen slice indices).

### Save semantics

- Prefer existing wrappers:
  - **Synthetic**: `createProjectStructure` → `saveDataset` → `saveRun`.
  - **Load**: `loadWorkspace` to load a saved `.mat` and automatically locate the matching `*_LOGfile.txt`.
- Minimum save payload for a run:
  - `log`, `data`, `params`, `meta` in a `.mat` (use `-v7.3` if large).
  - A copy of the active `*_LOGfile.txt` next to the `.mat`.

### UI / interactivity policy

- **Default mode is interactive** (intended for users).
- **Non-interactive mode is for reproduction only** (CI, paper figures, exact reruns).
- Scripts may be interactive, but interaction must be **localized**:
  - Prefer `uigetfile` / `uigetdir` inside wrappers (or one dedicated “Load/Select” block).
  - If `input(...)` is used, its outputs must be written into `params`/`cfg` and saved.
- For refactors, prefer adding `params.<stage>.interactive = true/false` (or `cfg.<stage>.interactive`) rather than scattering `input(...)` across the script.

#### UI principle

Interactive selection should feel like navigating a **checkpoint tree**:

- pick **project**
- pick **node/checkpoint**
- load it
- continue work → save as a **new node** (branch) with minimal friction

### Naming conventions

- Script names:
  - `run_<dataset>.m` for trunk orchestrators (thin).
  - Keep long, exploratory notebooks out of `run_*.m` (move into `examples/` or dedicated analysis scripts).
- Block IDs:
  - Follow the `run_synthetic_data.m` convention: `ABXXZ` (Category/Subcategory/Number/Variant).

### “Thin script, thick wrappers” rule

If a block grows beyond ~30–50 lines of logic (excluding plotting), create a wrapper function and call it from the script. This keeps refactors safe and makes blocks testable in isolation.

### Standard skeleton for new scripts (copy/paste)

```matlab
% TRUNK SCRIPT: <Title>
% =========================================================================
% Purpose: <one sentence>
%
%% SECTION 0: Path Initialization
clc; clear; close all;
repo_root = fileparts(fileparts(mfilename('fullpath')));
run(fullfile(repo_root, 'init_sbd.m'));

% Optional: initialize meta/log here or in the first block

%% GD01A: <Block summary>
% PRESETS
% params.<...> = ...;
% DO NOT EDIT BELOW
% ...

%% DONE
```

### Open questions for iteration (next design pass)

- For **real data**, should “project” live next to the raw data by default, or inside a centralized `projects/` root (like synthetic)?
- What is the minimal set of files that define each checkpoint node (so the UI can show a clean tree without scanning huge files)?

