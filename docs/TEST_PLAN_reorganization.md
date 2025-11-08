# MT-SBD-STM Reorganization Test Plan

**Date**: 2025-10-31  
**Author**: AI Assistant  
**Purpose**: Comprehensive testing of the three-phase workflow reorganization

---

## Overview

This document outlines a series of tests to validate that the reorganized MT-SBD-STM workflow achieves the intended goals:

1. **Hierarchical data organization** (`data.synGen`, `data.mcsbd_slice`, `data.mcsbd_block`)
2. **Automatic kernel initialization** with optional manual override
3. **Project-based folder structure** with organized saving
4. **Sequential loading** with automatic merging
5. **Three-phase workflow** with auto-save at boundaries

---

## Test Categories

### Category 1: Data Structure Validation
### Category 2: Kernel Initialization
### Category 3: File Saving and Organization
### Category 4: Sequential Loading
### Category 5: Workflow Integration
### Category 6: Edge Cases and Error Handling

---

## Test Category 1: Data Structure Validation

### Test 1.1: Synthetic Data Generation Structure
**Objective**: Verify that `generateSyntheticData` creates correct `data.synGen` structure

**Steps**:
1. Run `generateSyntheticData` with test parameters
2. Check that `data.synGen` exists and is a struct
3. Verify all expected fields exist: `Y`, `A0`, `A0_noiseless`, `X0`, `Y_ref`, `X0_ref`, `A0_ref`
4. Verify no `data.ground_truth` field exists (old structure)

**Expected Result**:
- `data.synGen` contains all synthetic generation outputs
- Field dimensions match parameters (e.g., `Y` is `[H×W×num_slices]`)
- No legacy field names present

**Pass Criteria**: All fields present with correct dimensions, no errors

---

### Test 1.2: Auto Kernel Initialization Structure
**Objective**: Verify that `autoInitializeKernels` creates correct `data.synGen.kernel_guess` structure

**Steps**:
1. Run `generateSyntheticData` followed by `autoInitializeKernels`
2. Check that `data.synGen.kernel_guess` exists
3. Verify fields: `kernel_centers`, `A_init`, `isolation_scores`, `method`
4. Verify `method` field equals `'auto'`

**Expected Result**:
- `kernel_guess.method` = `'auto'`
- `kernel_guess.A_init` is a cell array of initialized kernels
- `kernel_guess.kernel_centers` is `[K×2]` matrix

**Pass Criteria**: All fields present, method is 'auto', no errors

---

### Test 1.3: Manual Kernel Initialization Structure
**Objective**: Verify that `initializeKernelsRef` overwrites auto init correctly

**Steps**:
1. Run `generateSyntheticData` + `autoInitializeKernels` (creates auto init)
2. Run `initializeKernelsRef` (manual initialization)
3. Check that `data.synGen.kernel_guess.method` now starts with `'manual'`
4. Verify format is `'manual01'`
5. Run `initializeKernelsRef` again
6. Verify method increments to `'manual02'`

**Expected Result**:
- First manual init: `method` = `'manual01'`
- Second manual init: `method` = `'manual02'`
- Auto init is completely overwritten

**Pass Criteria**: Sequential numbering works, method field updates correctly

---

### Test 1.4: Algorithm Result Structure
**Objective**: Verify that algorithm results are stored in `data.mcsbd_slice`

**Steps**:
1. Run complete workflow through `decomposeReferenceSlice`
2. Check that `data.mcsbd_slice` exists
3. Verify fields: `A`, `X`, `b`, `extras`, `mtsbd_time`, `final_metrics`, `final_kernel_quality`
4. Verify no `data.singleslice_run` field exists (old structure)

**Expected Result**:
- `data.mcsbd_slice` contains all algorithm outputs
- Field dimensions match expected (e.g., `X` is `[H×W×K]`)
- No legacy field names

**Pass Criteria**: All fields present with correct dimensions, no errors

---

## Test Category 2: Kernel Initialization

### Test 2.1: Auto Initialization from Ground Truth
**Objective**: Verify auto initialization works without running MT-SBD

**Steps**:
1. Run `generateSyntheticData`
2. Run `autoInitializeKernels`
3. Verify that kernels are initialized at isolated points
4. Check that no MT-SBD algorithm was run (no `data.mcsbd_slice`)

**Expected Result**:
- Kernels initialized successfully
- Isolation scores calculated
- No algorithm results present yet

**Pass Criteria**: Initialization completes, kernels look reasonable, no algorithm run

---

### Test 2.2: Manual Initialization Override
**Objective**: Verify manual initialization properly overwrites auto

**Steps**:
1. Run auto initialization (creates `auto` method)
2. Note the `A_init` kernels
3. Run manual initialization
4. Verify `A_init` kernels are different
5. Verify method changed to `manual01`

**Expected Result**:
- Manual kernels replace auto kernels
- Method field updates correctly
- Previous auto init is completely overwritten

**Pass Criteria**: Kernels are different, method updates, no remnants of auto init

---

### Test 2.3: Multiple Manual Iterations
**Objective**: Verify sequential numbering for multiple manual initializations

**Steps**:
1. Run auto init
2. Run manual init → should be `manual01`
3. Run manual init again → should be `manual02`
4. Run manual init again → should be `manual03`

**Expected Result**:
- Each manual init increments the counter
- Method follows pattern: `manual01`, `manual02`, `manual03`, etc.

**Pass Criteria**: Sequential numbering works correctly, no skips or duplicates

---

## Test Category 3: File Saving and Organization

### Test 3.1: Project Structure Creation
**Objective**: Verify `createProjectStructure` creates correct folder hierarchy

**Steps**:
1. Call `createProjectStructure()`
2. Check that `projects/` folder is created
3. Check that `synthetic_<timestamp>/` subfolder is created
4. Verify `meta` struct contains correct paths

**Expected Result**:
- Folder structure: `\projects\synthetic_YYYYMMDD_HHMMSS\`
- `meta.project_path` points to created folder
- `meta.timestamp` matches folder name

**Pass Criteria**: Folders created, meta struct populated correctly

---

### Test 3.2: Auto Dataset Saving
**Objective**: Verify `saveDataset` saves auto-initialized dataset correctly

**Steps**:
1. Run workflow through auto initialization
2. Call `saveDataset(log, data, params, meta)`
3. Check that `auto.mat` file is created in project folder
4. Check that `auto_LOGfile.txt` is created
5. Load `auto.mat` and verify contents

**Expected Result**:
- Files created: `auto.mat`, `auto_LOGfile.txt`
- `auto.mat` contains `log`, `data`, `params`
- `data.synGen.kernel_guess.method` = `'auto'`

**Pass Criteria**: Files created, contents correct, no errors

---

### Test 3.3: Manual Dataset Saving
**Objective**: Verify `saveDataset` saves manual-initialized dataset with correct naming

**Steps**:
1. Run workflow through manual initialization (first time)
2. Call `saveDataset(log, data, params, meta)`
3. Check that `manual01.mat` is created
4. Run manual initialization again (second time)
5. Call `saveDataset` again
6. Check that `manual02.mat` is created

**Expected Result**:
- First save: `manual01.mat`, `manual01_LOGfile.txt`
- Second save: `manual02.mat`, `manual02_LOGfile.txt`
- Both files contain correct `method` field

**Pass Criteria**: Sequential numbering works, files created correctly

---

### Test 3.4: Run Results Saving
**Objective**: Verify `saveRun` creates runXX subfolder structure

**Steps**:
1. Complete workflow through algorithm execution
2. Call `saveRun(log, data, params, meta)`
3. Check that `run01/` folder is created
4. Check that `run01/run01.mat` and `run01/run01_LOGfile.txt` are created
5. Run algorithm again with different parameters
6. Call `saveRun` again
7. Check that `run02/` folder is created

**Expected Result**:
- First save: `run01/run01.mat`, `run01/run01_LOGfile.txt`
- Second save: `run02/run02.mat`, `run02/run02_LOGfile.txt`
- Files contain algorithm results in `data.mcsbd_slice`

**Pass Criteria**: Folders and files created, sequential numbering works

---

### Test 3.5: Folder Structure Integrity
**Objective**: Verify complete folder structure after full workflow

**Steps**:
1. Run complete workflow with auto init
2. Run manual init and save
3. Run algorithm and save
4. Check folder structure

**Expected Structure**:
```
\projects\synthetic_YYYYMMDD_HHMMSS\
  ├─ auto.mat
  ├─ auto_LOGfile.txt
  ├─ manual01.mat
  ├─ manual01_LOGfile.txt
  └─ run01\
     ├─ run01.mat
     └─ run01_LOGfile.txt
```

**Pass Criteria**: All files and folders present in correct hierarchy

---

## Test Category 4: Sequential Loading

### Test 4.1: Load Dataset File
**Objective**: Verify `loadSequential` loads dataset file correctly

**Steps**:
1. Save a dataset file (`auto.mat`)
2. Clear workspace
3. Call `loadSequential('workspace_file', 'path/to/auto.mat')`
4. Verify `data.synGen` is loaded
5. Verify no `data.mcsbd_slice` exists

**Expected Result**:
- `data.synGen` loaded with all fields
- `params` and `log` loaded
- `meta.file_type` = `'dataset'`

**Pass Criteria**: Dataset loaded correctly, no algorithm results present

---

### Test 4.2: Load Run File with Auto-Merge
**Objective**: Verify `loadSequential` auto-loads parent dataset when loading run file

**Steps**:
1. Save dataset file (`auto.mat`)
2. Save run file (`run01/run01.mat`)
3. Clear workspace
4. Call `loadSequential('workspace_file', 'path/to/run01/run01.mat')`
5. Verify both `data.synGen` AND `data.mcsbd_slice` are present

**Expected Result**:
- `data.synGen` loaded from parent dataset
- `data.mcsbd_slice` loaded from run file
- Both substructs merged into single `data`
- `meta.load_history` shows both files loaded

**Pass Criteria**: Both data substructs present, merge successful

---

### Test 4.3: Load with Memory Protection
**Objective**: Verify memory overwrite protection works

**Steps**:
1. Load a workspace (creates `log`, `data`, `params`)
2. Call `loadSequential` again
3. Verify prompt appears asking to delete/rename/cancel
4. Choose "rename" with suffix "_old"
5. Verify new variables created: `log_old`, `data_old`, `params_old`
6. Verify new workspace loaded into `log`, `data`, `params`

**Expected Result**:
- Prompt appears correctly
- Old workspace preserved with suffix
- New workspace loaded into base variables

**Pass Criteria**: Memory protection works, no data loss

---

### Test 4.4: Load History Tracking
**Objective**: Verify `meta.load_history` tracks all loaded files

**Steps**:
1. Load run file (which auto-loads dataset)
2. Check `meta.load_history.files_loaded`
3. Check `meta.load_history.load_order`

**Expected Result**:
- `files_loaded` contains paths to both dataset and run files
- `load_order` shows dataset loaded first, then run file
- Correct chronological order

**Pass Criteria**: Load history accurate and complete

---

## Test Category 5: Workflow Integration

### Test 5.1: Complete Auto Workflow
**Objective**: Run complete workflow with auto initialization only

**Steps**:
1. Run `run_synthetic_data.m` with `run_manual_init = false`
2. Verify auto initialization runs
3. Verify Phase 1 save creates `auto.mat`
4. Verify algorithm blocks run
5. Verify Phase 2 save creates `run01/run01.mat`
6. Check all data structures are correct

**Expected Result**:
- Workflow completes without errors
- Files created: `auto.mat`, `run01/run01.mat`
- Data structures correct at each phase

**Pass Criteria**: Complete workflow runs, all saves occur, data correct

---

### Test 5.2: Complete Manual Workflow
**Objective**: Run complete workflow with manual initialization

**Steps**:
1. Run `run_synthetic_data.m` with `run_manual_init = true`
2. Perform manual kernel selection
3. Verify Phase 1 save creates `manual01.mat`
4. Verify algorithm blocks run
5. Verify Phase 2 save creates `run01/run01.mat`

**Expected Result**:
- Manual initialization prompts appear
- Files created: `manual01.mat`, `run01/run01.mat`
- Manual kernels used in algorithm

**Pass Criteria**: Manual init works, saves correct, algorithm uses manual kernels

---

### Test 5.3: Multiple Run Workflow
**Objective**: Run algorithm multiple times with different parameters

**Steps**:
1. Complete workflow once (creates `run01`)
2. Modify algorithm parameters
3. Run Phase 2 blocks again
4. Verify `run02/` folder created
5. Repeat with different parameters
6. Verify `run03/` folder created

**Expected Result**:
- Sequential run folders: `run01/`, `run02/`, `run03/`
- Each contains independent results
- No overwriting of previous runs

**Pass Criteria**: Multiple runs saved independently, sequential numbering works

---

### Test 5.4: Resume from Saved Dataset
**Objective**: Load saved dataset and continue to algorithm phase

**Steps**:
1. Run Phase 1 only (generate + init + save)
2. Close MATLAB / clear workspace
3. Load saved dataset using `loadSequential`
4. Run Phase 2 blocks (algorithm execution)
5. Verify results saved correctly

**Expected Result**:
- Dataset loads correctly
- Algorithm runs on loaded data
- Run results saved in correct location

**Pass Criteria**: Resume works, no data loss, saves in correct project folder

---

### Test 5.5: Resume from Saved Run
**Objective**: Load saved run and continue to visualization

**Steps**:
1. Complete Phase 1 and Phase 2 (save run)
2. Close MATLAB / clear workspace
3. Load saved run using `loadSequential`
4. Verify both `data.synGen` and `data.mcsbd_slice` present
5. Run Phase 3 (visualization)

**Expected Result**:
- Both data substructs loaded
- Visualization works with loaded data
- No errors accessing any fields

**Pass Criteria**: Resume works, all data accessible, visualization completes

---

## Test Category 6: Edge Cases and Error Handling

### Test 6.1: Missing Parent Dataset
**Objective**: Test behavior when run file's parent dataset is missing

**Steps**:
1. Create and save run file
2. Delete parent dataset file
3. Attempt to load run file with `loadSequential`

**Expected Result**:
- Warning message about missing dataset
- Run file still loads (with only `data.mcsbd_slice`)
- No crash

**Pass Criteria**: Graceful handling, warning issued, partial load succeeds

---

### Test 6.2: Corrupted File Loading
**Objective**: Test behavior when loading corrupted/invalid file

**Steps**:
1. Create invalid .mat file (wrong structure)
2. Attempt to load with `loadSequential`

**Expected Result**:
- Clear error message about invalid file
- No workspace corruption
- User can retry with different file

**Pass Criteria**: Error caught, clear message, no crash

---

### Test 6.3: Disk Space Issues
**Objective**: Test behavior when disk is full during save

**Steps**:
1. (Simulated) Fill disk to near capacity
2. Attempt to save large workspace

**Expected Result**:
- Error message about disk space
- Partial file cleaned up
- User notified to free space

**Pass Criteria**: Error caught, no partial files left, clear message

---

### Test 6.4: Concurrent Access
**Objective**: Test behavior when multiple MATLAB instances access same project

**Steps**:
1. Open project in two MATLAB instances
2. Attempt to save from both simultaneously

**Expected Result**:
- File locking or sequential numbering prevents collision
- Both saves succeed with different run numbers
- No data corruption

**Pass Criteria**: No collisions, both saves succeed, data intact

---

## Test Execution Checklist

### Phase 1: Basic Functionality
- [ ] Test 1.1: Data structure (synGen)
- [ ] Test 1.2: Auto init structure
- [ ] Test 1.3: Manual init structure
- [ ] Test 1.4: Algorithm result structure
- [ ] Test 2.1: Auto initialization
- [ ] Test 2.2: Manual override
- [ ] Test 3.1: Project structure creation
- [ ] Test 3.2: Auto dataset saving
- [ ] Test 3.3: Manual dataset saving

### Phase 2: Advanced Functionality
- [ ] Test 2.3: Multiple manual iterations
- [ ] Test 3.4: Run results saving
- [ ] Test 3.5: Folder structure integrity
- [ ] Test 4.1: Load dataset file
- [ ] Test 4.2: Load run file with merge
- [ ] Test 4.3: Memory protection
- [ ] Test 4.4: Load history tracking

### Phase 3: Integration Testing
- [ ] Test 5.1: Complete auto workflow
- [ ] Test 5.2: Complete manual workflow
- [ ] Test 5.3: Multiple run workflow
- [ ] Test 5.4: Resume from dataset
- [ ] Test 5.5: Resume from run

### Phase 4: Edge Cases
- [ ] Test 6.1: Missing parent dataset
- [ ] Test 6.2: Corrupted file loading
- [ ] Test 6.3: Disk space issues
- [ ] Test 6.4: Concurrent access

---

## Success Criteria

The reorganization is considered successful if:

1. **All Category 1-2 tests pass** (Data structure and initialization)
2. **All Category 3 tests pass** (File saving and organization)
3. **All Category 4 tests pass** (Sequential loading)
4. **At least 4/5 Category 5 tests pass** (Workflow integration)
5. **At least 3/4 Category 6 tests pass** (Edge cases)

---

## Notes for Test Execution

1. **Test Environment**: Use a separate test directory to avoid affecting production data
2. **Test Data**: Use small synthetic datasets (N_obs=20, num_slices=2) for speed
3. **Logging**: Keep detailed logs of each test execution
4. **Regression**: After fixing any issues, re-run all previous tests
5. **Documentation**: Update this plan as new edge cases are discovered

---

**End of Test Plan**

