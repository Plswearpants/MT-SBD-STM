# Progress Summary - Session 2025-11-04

## 🎯 Main Achievement: Unified Architecture for Params Management

### **Core Principle Established:**
**"Flat internally, hierarchical for storage"**
- All wrapper functions work with **flat** `params` internally
- All saved files store **hierarchical** `params`
- The script **NEVER** accesses `params` fields directly
- Conversion handled by `organizeParams.m` utility

---

## ✅ Completed Work

### 1. **Created `organizeParams.m` Utility** ✅
- **Location**: `Dong_func/wrapper/organizeParams.m`
- **Purpose**: Bidirectional conversion between flat and hierarchical params
- **Modes**:
  - `'write'`: flat → hierarchical (for storage)
  - `'extract'`: hierarchical → flat (for processing)
- **Field Mapping**:
  - `params.synGen.*` (e.g., SNR, num_kernels, ref_slice)
  - `params.mcsbd_slice.*` (e.g., lambda1, maxIT)
  - `params.mcsbd_block.*` (future use)

### 2. **Updated Wrapper Functions** ✅

#### ✅ `generateSyntheticData.m`
- Added `log` as first required parameter
- Uses flat params internally
- Converts to hierarchical at end with `organizeParams(params, 'write')`
- **Moved logging inside function**:
  - Block start: "GD01A"
  - Generation details logging
- Updated documentation

#### ✅ `autoInitializeKernels.m`
- Added `log` as first required parameter
- Extracts to flat at start (if hierarchical input)
- Uses flat params internally
- Converts to hierarchical at end
- **Moved logging inside function**:
  - Auto-initialization results logging
- Updated documentation

#### ✅ `initializeKernelsRef.m`
- Added `log` as first required parameter
- Extracts to flat at start (if hierarchical input)
- Uses flat params internally
- Converts to hierarchical at end
- **Moved logging inside function**:
  - Manual initialization results logging
- Updated documentation

### 3. **Created `autoSave.m` Wrapper** ✅
- **Location**: `Dong_func/wrapper/autoSave.m`
- **Purpose**: Intelligent automatic saving based on workflow phase
- **Features**:
  - Auto-detects phase from `data` struct:
    - Pre-run: `data.synGen.kernel_guess` exists → calls `saveDataset`
    - Post-run: `data.mcsbd_slice` or `data.mcsbd_block` exists → calls `saveRun`
  - Can force specific phase with `'phase'` parameter
  - Verbose output control
- **Note**: Currently NOT used in script (reverted to direct `saveDataset`/`saveRun` calls per user preference)

### 4. **Updated `run_synthetic_data.m` Script** ✅

#### **Removed ALL `organizeParams(..., 'extract')` calls from script**
- Script never accesses `params` fields directly
- All field access happens inside wrapper functions
- Cleaner, more maintainable architecture

#### **Kept `organizeParams(..., 'write')` before saves** (3 locations)
1. **Line 121**: GD01A block - before auto initialization save
2. **Line 172**: IK01A block - before manual initialization save
3. **Line 538**: Phase 2 completion - before run save

#### **Updated Block Structure** (3 completed blocks)
New pattern for all blocks:
```matlab
% -------------------------------------------------------------------------
% PRESETS: User-configurable parameters
% -------------------------------------------------------------------------
param1 = value1;
param2 = value2;

%%%%%%%%%%%%%%%%%% DO NOT EDIT BELOW %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Function call (wrapper handles its own logging)
[data, params] = wrapperFunction(log, data, params, ...
    'param1', param1, 'param2', param2);

% Clear preset variables
clearvars param1 param2
```

**Blocks Updated:**
- ✅ GD01A: Generate data + auto initialization
- ✅ IK01A: Manual initialization (optional)
- ✅ LD01A: Load dataset (removed extract after load)

### 5. **Fixed Logging Path Bug** 🐛
- **Issue**: `logUsedBlocks.m` used `strcat(LOGpath, '/', LOGfile, ...)` which failed on Windows
- **Note**: Fix was reverted by user (will need to reapply if issue persists)
- **Solution**: Use `fullfile(LOGpath, [LOGfile '_LOGfile.txt'])` for cross-platform compatibility

### 6. **Fixed log.file Cell Array Bug** 🐛

#### Root Cause Analysis:
- **Issue**: `log.file` was sometimes a **string object** instead of **char array**, causing type inconsistencies
- **When it happened**:
  - ✅ User types custom name → `input(prompt, "s")` → **char array** (correct)
  - ❌ User presses Enter (default) → `defaultName = "Project"` → **string object** (incorrect)
  - ❌ Unique name appended → `strcat(nameString, "_", sprintf("%03d", n))` → **string object** (incorrect)

#### Symptoms:
- **In save functions**: `char(fullfile(...))` created 2D char array:
  ```matlab
  ans(:,:,1) = 'C:\Users\...\MT-SBD-STM\tests\NewProject_005'
  ans(:,:,2) = 'C:\Users\...\MT-SBD-STM\tests\_LOGfile.txt  '
  ```

#### Solutions Applied:

**Solution 1 (Root Fix)**: Modified `setLogFile.m` to always return char
```matlab
% After user input or default assignment
nameString = char(nameString);  % Force conversion to char

% In uniqueness loop
testNameString = char(strcat(...));  % Ensure appended names are char
nameString = char(testNameString);   // Final assignment as char
```
- **Status**: ✅ Fixed in `Dong_func/basic/setLogFile.m`

**Solution 2 (Defensive Fix)**: Added type handling in save functions
```matlab
% Handle log.file (convert cell/string to char if necessary)
if iscell(log.file)
    log_file_str = char(log.file{1});
else
    log_file_str = char(log.file);
end
log_source = fullfile(log.path, [log_file_str '_LOGfile.txt']);
```
- **Status**: ✅ Fixed in `saveDataset.m` and `saveRun.m` (defensive programming)

**Why both fixes?**:
- Root fix ensures `log.file` is always char from the start
- Defensive fix protects against edge cases when loading old workspaces with string-type `log.file`

---

## 📋 TODO List for Tomorrow

### **Priority 1: Complete Remaining Wrapper Functions** 🔴

#### 1. **Update `decomposeReferenceSlice.m`** ⏳
- [ ] Add `log` as first required parameter
- [ ] Already has extract at start (line 56-58) ✅
- [ ] Already has write at end (line 236) ✅
- [ ] **Move logging inside function**:
  - [ ] Block start logging ("DS01A")
  - [ ] Phase I/II logging
  - [ ] Convergence logging
- [ ] Update corresponding **DS01A block** in script:
  - [ ] Add `log` parameter to function call
  - [ ] Remove external logging
  - [ ] Verify block structure (presets → function → clearvars)

#### 2. **Update `findIsolatedPoints.m`** ⏳
- [ ] Add `log` as first required parameter
- [ ] Already has extract at start (line 69-71) ✅
- [ ] **Add write at end** ❌ (MISSING!)
- [ ] **Move logging inside function**:
  - [ ] Isolation analysis results
- [ ] Update corresponding **IS01A block** in script:
  - [ ] Add `log` parameter to function call
  - [ ] Remove external logging
  - [ ] Verify block structure

#### 3. **Update Other Wrapper Functions** (if any used in script)
- [ ] Check which other wrappers are called in script
- [ ] Apply same pattern to each
- [ ] Verify logging is moved inside

### **Priority 2: Testing & Validation** 🟡

#### 1. **Test Full Workflow** 🧪
- [ ] **Test Phase 1 (Pre-run)**:
  - [ ] Run GD01A: generate data + auto init + save
  - [ ] Verify `auto.mat` is created with hierarchical params
  - [ ] Run IK01A (optional): manual init + save
  - [ ] Verify `manual01.mat` is created
  - [ ] Check log files are copied correctly

- [ ] **Test Phase 2 (Post-run)**:
  - [ ] Run LD01A: load dataset
  - [ ] Verify params are loaded as hierarchical
  - [ ] Run DS01A: decompose reference slice
  - [ ] Run IS01A: find isolated points
  - [ ] Verify run results saved correctly

- [ ] **Test Param Conversion**:
  - [ ] Load saved workspace
  - [ ] Verify params structure is hierarchical
  - [ ] Run wrapper function
  - [ ] Verify wrapper can extract to flat
  - [ ] Verify wrapper returns hierarchical

#### 2. **Edge Case Testing**
- [ ] Test loading workspace with no prior workspace
- [ ] Test running IK01A multiple times (manual01, manual02, ...)
- [ ] Test running algorithm on different datasets (auto vs manual)
- [ ] Test memory protection when loading

### **Priority 3: Documentation & Cleanup** 🟢

#### 1. **Update Documentation**
- [ ] Update README with new architecture
- [ ] Document `organizeParams` usage patterns
- [ ] Add examples of wrapper function patterns
- [ ] Update folder structure documentation

#### 2. **Code Review & Cleanup**
- [ ] Check for any orphaned `organizeParams` calls
- [ ] Verify all wrappers follow consistent pattern
- [ ] Review all logging calls (should be inside wrappers)
- [ ] Clean up any temporary files/comments

#### 3. **Consider Future Enhancements**
- [ ] Decide whether to use `autoSave` or keep explicit `saveDataset`/`saveRun`
- [ ] Review if `autoSave` phase detection logic is robust enough
- [ ] Consider adding validation to `organizeParams` (check for missing fields)

---

## 🏗️ Architecture Decisions Made

### **1. Params Conversion Pattern**

**In Wrapper Functions:**
```matlab
function [data, params] = wrapperFunction(log, data, params, ...)
    % Parse inputs
    parse(p, log, data, params, varargin{:});
    
    % Extract to flat (if hierarchical input)
    if isfield(params, 'synGen')
        params = organizeParams(params, 'extract');
    end
    
    % Work with flat params internally
    some_value = params.SNR;  % Direct access OK inside function
    params.new_field = 123;    % Direct modification OK
    
    % LOG: Internal logging
    LOGcomment = sprintf("...", params.num_kernels);
    LOGcomment = logUsedBlocks(log.path, log.file, "BLOCKID", LOGcomment, 0);
    
    % Convert to hierarchical for output
    params = organizeParams(params, 'write');
end
```

**In Script:**
```matlab
% NEVER access params fields directly in script!
% Only pass params to functions and receive back

% Before saving: ensure hierarchical
params = organizeParams(params, 'write');
meta = saveDataset(log, data, params, meta);
% NO extraction after save!

% After loading: params is already hierarchical
[log, data, params, meta] = loadWorkspace(...);
% NO extraction needed - next wrapper will extract internally
```

### **2. Logging Pattern**

**OLD (in script):**
```matlab
LOGcomment = sprintf("Details");
LOGcomment = logUsedBlocks(log.path, log.file, "BLOCKID", LOGcomment, 0);
[data, params] = wrapperFunction(...);
LOGcomment = sprintf("Results");
LOGcomment = logUsedBlocks(log.path, log.file, "  ^  ", LOGcomment, 0);
```

**NEW (in wrapper):**
```matlab
% In script - only log presets
LOGcomment = sprintf("Presets: param=%s", param);
LOGcomment = logUsedBlocks(log.path, log.file, "BLOCKID", LOGcomment, 0);

% Function handles its own result logging
[data, params] = wrapperFunction(log, data, params, ...);
```

### **3. Block Structure in Script**

**Standard Pattern:**
```matlab
%% =========================================================================
%% BLOCKID: Block-Name-Description
%  =========================================================================
%  Description of what the block does
%
%  Prerequisites: ...
%  Dependencies: ...

% -------------------------------------------------------------------------
% PRESETS: User-configurable parameters
% -------------------------------------------------------------------------
param1 = value1;
param2 = value2;

%%%%%%%%%%%%%%%%%% DO NOT EDIT BELOW %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% LOG: Block start (preset parameters only)
LOGcomment = sprintf("Presets: param1=%s, param2=%s", param1, param2);
LOGcomment = logUsedBlocks(log.path, log.file, "BLOCKID", LOGcomment, 0);

% Function call (compact, one-liner preferred, wrapper handles internal logging)
[data, params] = wrapperFunction(log, data, params, ...
    'param1', param1, 'param2', param2);

% Clear preset variables
clearvars param1 param2
```

---

## 📊 Current Status

### **Wrapper Functions Status:**
| Function | Log Added | Extract Start | Write End | Script Updated | Status |
|----------|-----------|---------------|-----------|----------------|--------|
| `generateSyntheticData` | ✅ | N/A (generator) | ✅ | ✅ | **DONE** |
| `autoInitializeKernels` | ✅ | ✅ | ✅ | ✅ | **DONE** |
| `initializeKernelsRef` | ✅ | ✅ | ✅ | ✅ | **DONE** |
| `decomposeReferenceSlice` | ✅ | ✅ | ✅ | ✅ | **DONE** |
| `findIsolatedPoints` | ✅ | ✅ | ✅ | ✅ | **DONE** |

### **Files Modified Today:**
1. ✅ `Dong_func/wrapper/organizeParams.m` - **CREATED**
2. ✅ `Dong_func/wrapper/autoSave.m` - **CREATED**
3. ✅ `Dong_func/wrapper/generateSyntheticData.m` - **UPDATED** (log param, internal logging)
4. ✅ `Dong_func/wrapper/autoInitializeKernels.m` - **UPDATED** (log param, internal logging)
5. ✅ `Dong_func/wrapper/initializeKernelsRef.m` - **UPDATED** (log param, internal logging)
6. ✅ `Dong_func/wrapper/decomposeReferenceSlice.m` - **UPDATED** (log param, internal logging)
7. ✅ `Dong_func/wrapper/findIsolatedPoints.m` - **UPDATED** (log param, internal logging, write at end, ref_results input)
8. ✅ `Dong_func/wrapper/saveDataset.m` - **UPDATED** (defensive fix for log.file)
9. ✅ `Dong_func/wrapper/saveRun.m` - **UPDATED** (defensive fix for log.file)
10. ✅ `Dong_func/basic/setLogFile.m` - **UPDATED** (root fix for log.file string/char issue)
11. ✅ `scripts/run_synthetic_data.m` - **UPDATED** (GD01A, IK01A, LD01A, DS01A, IS01A blocks)
12. 🔄 `Dong_func/basic/logUsedBlocks.m` - **ATTEMPTED FIX** (reverted by user)

---

## 🚀 Quick Start for Tomorrow

### **Step 1: Update `decomposeReferenceSlice.m`**
```bash
# Open file
code Dong_func/wrapper/decomposeReferenceSlice.m

# Follow pattern from generateSyntheticData.m:
# 1. Add 'log' parameter
# 2. Move logging inside
# 3. Verify extract at start
# 4. Verify write at end
```

### **Step 2: Update DS01A Block in Script**
```bash
# Open file
code scripts/run_synthetic_data.m

# Find DS01A block
# Update function call to include 'log'
# Remove external logging
```

### **Step 3: Repeat for `findIsolatedPoints.m`**
Same pattern as Step 1 & 2

### **Step 4: Test**
Run the full workflow in MATLAB and verify:
- Saves work correctly
- Loads work correctly
- Logging works correctly
- No errors!

---

## 📝 Notes & Reminders

1. **Windows Path Issue**: If log file warnings appear, reapply the `fullfile()` fix to `logUsedBlocks.m`

2. **autoSave Usage**: Currently NOT using `autoSave.m` in script - keeping explicit `saveDataset`/`saveRun` calls per user preference. Decision can be revisited later.

3. **Critical Architecture Rule**: 
   - **Script = orchestration only** (no direct params access)
   - **Wrappers = processing + logging** (handle their own conversions)
   - **Storage = always hierarchical** (consistent file format)

4. **Logging Strategy**:
   - Block start (in script): Log preset parameters
   - Block results (in wrapper): Log execution results
   - This keeps script clean while maintaining complete logs

---

## 🎓 Key Lessons Learned

1. **"Flat internally, hierarchical for storage"** creates clean separation
2. **Never extract after save** - script doesn't need params after that point
3. **Wrappers should be self-contained** - handle their own logging and conversions
4. **Use `fullfile()` not `strcat()` for paths** - cross-platform compatibility
5. **Script orchestrates, wrappers process** - clear separation of concerns

---

**End of Session 2025-11-04**

