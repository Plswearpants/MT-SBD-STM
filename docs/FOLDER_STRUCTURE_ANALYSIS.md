# Folder Structure Analysis: Current vs. Proposed

## Current Structure

```
projects/
  synthetic_<timestamp>/
    ├─ auto.mat                          (dataset file)
    ├─ auto_LOGfile.txt                  (dataset log)
    ├─ manual01.mat                      (dataset file)
    ├─ manual01_LOGfile.txt              (dataset log)
    ├─ manual02.mat                      (dataset file)
    ├─ manual02_LOGfile.txt              (dataset log)
    ├─ run_auto_01/                      (run folder)
    │  ├─ run_auto_01.mat
    │  └─ run_auto_01_LOGfile.txt
    ├─ run_auto_02/
    │  ├─ run_auto_02.mat
    │  └─ run_auto_02_LOGfile.txt
    ├─ run_manual01_01/
    │  ├─ run_manual01_01.mat
    │  └─ run_manual01_01_LOGfile.txt
    └─ run_manual02_01/
       ├─ run_manual02_01.mat
       └─ run_manual02_01_LOGfile.txt
```

**Characteristics:**
- Flat structure for dataset files (all .mat files at project root)
- Run folders at project root with dataset name in folder name
- Visualization files saved in `log.path` (typically `tests/` directory), not organized in project structure
- Dataset name embedded in run folder name: `run_<dataset>_XX`

## Proposed Nested Structure

```
projects/
  synthetic_<timestamp>/
    ├─ auto/                              (dataset folder)
    │  ├─ auto.mat                        (dataset file)
    │  ├─ auto_LOGfile.txt                (dataset log)
    │  ├─ run01/                          (run folder)
    │  │  ├─ run01.mat
    │  │  ├─ run01_LOGfile.txt
    │  │  ├─ visualization01/             (visualization folder)
    │  │  │  ├─ figure1.fig
    │  │  │  ├─ figure2.png
    │  │  │  └─ ...
    │  │  ├─ visualization02/
    │  │  │  └─ ...
    │  │  └─ visualization03/
    │  │     └─ ...
    │  ├─ run02/
    │  │  ├─ run02.mat
    │  │  ├─ run02_LOGfile.txt
    │  │  └─ visualization01/
    │  │     └─ ...
    │  └─ run03/
    │     └─ ...
    ├─ manual01/                          (dataset folder)
    │  ├─ manual01.mat
    │  ├─ manual01_LOGfile.txt
    │  ├─ run01/
    │  │  ├─ run01.mat
    │  │  ├─ run01_LOGfile.txt
    │  │  └─ visualization01/
    │  │     └─ ...
    │  └─ run02/
    │     └─ ...
    └─ manual02/
       ├─ manual02.mat
       ├─ manual02_LOGfile.txt
       └─ run01/
          └─ ...
```

**Characteristics:**
- Hierarchical nesting: dataset → run → visualization
- Each level clearly indicates workflow phase
- Simple run numbering (run01, run02) since dataset is parent folder
- Visualization results organized within each run

---

## Pros and Cons Analysis

### Current Structure

#### Pros:
1. **Flat dataset access**: Easy to see all datasets at a glance
2. **Simple loading**: Direct path to dataset file (e.g., `auto.mat`)
3. **Clear run naming**: Dataset name embedded in run folder (`run_auto_01`)
4. **Less nesting**: Fewer folder levels to navigate
5. **Backward compatible**: Existing code expects files at project root

#### Cons:
1. **Cluttered project root**: Mix of files and folders at same level
2. **No visualization organization**: Visualization files scattered in separate directory
3. **Run naming redundancy**: Dataset name repeated in every run folder name
4. **Harder to browse dataset variants**: All dataset files mixed together
5. **No logical grouping**: Related files (dataset + runs) not grouped together

---

### Proposed Nested Structure

#### Pros:
1. **Clear workflow hierarchy**: Dataset → Run → Visualization reflects logical order
2. **Better organization**: Related files grouped together (dataset + all its runs)
3. **Cleaner project root**: Only dataset folders visible at top level
4. **Easier dataset browsing**: Each dataset folder contains all related work
5. **Simpler run numbering**: Just `run01`, `run02` (no need for dataset prefix)
6. **Visualization organization**: Visualizations properly nested within runs
7. **Scalability**: Easy to add new organizational levels (e.g., parameter sweeps)
8. **Self-documenting**: Folder structure tells the story of the workflow

#### Cons:
1. **More nesting**: Deeper folder hierarchy (dataset/run/visualization)
2. **Longer paths**: More characters in file paths
3. **Breaking change**: Requires updating all save/load functions
4. **Migration needed**: Existing projects would need restructuring
5. **More complex path logic**: Need to track dataset folder, run folder, visualization folder
6. **Loading complexity**: Need to navigate to dataset folder first

---

## Recommendation

**The proposed nested structure is better for the following reasons:**

1. **Better Scalability**: As projects grow, the nested structure prevents root-level clutter
2. **Logical Organization**: Matches the workflow: generate dataset → run algorithm → visualize
3. **Easier Dataset Management**: All work on a dataset is in one place
4. **Proper Visualization Organization**: Currently missing feature that nested structure provides
5. **Future-Proof**: Easy to extend (e.g., add parameter sweep folders, comparison folders)

**However, implementation considerations:**

1. **Migration Strategy**: Need to provide migration script for existing projects
2. **Path Updates**: Update all save functions (`saveDataset`, `saveRun`) and load functions
3. **Meta Struct Updates**: Track dataset folder path, run folder path, visualization folder path
4. **Backward Compatibility**: Consider compatibility mode for loading old structures

---

## Implementation Impact

### Functions to Update:
1. **`createProjectStructure.m`**: Create dataset folders instead of files at root
2. **`saveDataset.m`**: Save to `meta.project_path/<dataset>/<dataset>.mat`
3. **`saveRun.m`**: Save to `meta.project_path/<dataset>/runXX/runXX.mat`
4. **`loadWorkspace.m`**: Navigate nested structure (already uses UI, so minimal change)
5. **`autoSave.m`**: Update path logic for nested structure
6. **Visualization saving**: Create `visualizationXX/` folders within runs

### Meta Struct Changes:
```matlab
meta.project_path           % projects/synthetic_<timestamp>/
meta.dataset_folder         % auto/ or manual01/
meta.dataset_path           % projects/.../auto/auto.mat
meta.run_folder            % run01/
meta.run_path              % projects/.../auto/run01/
meta.visualization_folder  % visualization01/
meta.visualization_path     % projects/.../auto/run01/visualization01/
```

---

## Decision Matrix

| Criteria | Current | Proposed | Winner |
|----------|---------|----------|--------|
| Organization | ⭐⭐ | ⭐⭐⭐⭐⭐ | Proposed |
| Ease of Loading | ⭐⭐⭐⭐ | ⭐⭐⭐ | Current |
| Scalability | ⭐⭐ | ⭐⭐⭐⭐⭐ | Proposed |
| Visualization Organization | ⭐ | ⭐⭐⭐⭐⭐ | Proposed |
| Implementation Complexity | ⭐⭐⭐⭐ | ⭐⭐ | Current |
| Self-Documentation | ⭐⭐ | ⭐⭐⭐⭐⭐ | Proposed |
| Backward Compatibility | ⭐⭐⭐⭐⭐ | ⭐ | Current |

**Overall: Proposed structure wins for long-term maintainability and organization**

