# Template Documentation: Modular Data Processing Structure

## Overview

This template provides a standardized structure for organizing complex data processing workflows with comprehensive logging and reproducibility. It is based on the UBC_LAIR mainTrunk.m architecture.

---

## Core Principles

### 1. **Block-Based Architecture**
- Each functional unit is a self-contained "block"
- Blocks are organized by category (Loading, Processing, Selecting, Visualizing, etc.)
- Each block has a unique identifier for tracking and logging

### 2. **Comprehensive Logging**
- Every operation is logged with timestamp and parameters
- Logs enable full reproducibility of analysis
- Saved outputs are linked to their corresponding log entries

### 3. **Standardized Structure**
- Consistent format across all blocks
- Clear separation of user presets and execution code
- Automatic cleanup of temporary variables

---

## Block Structure

### Block Identifier System: `ABXXZ`

#### Component Breakdown:

**A - Category (1st character):**
- `L` = Loading (reading/importing data)
- `P` = Processing (data transformation, filtering)
- `V` = Visualizing (plotting/generating figures)
- `S` = Selecting (masking, region selection)
- `W` = Writing/Saving (exporting data, saving workspace)

**B - Subcategory (2nd character):**
- `A` = Averaging
- `D` = Data/Derivative
- `F` = Flatten
- `G` = Grid
- `I` = Initialize
- `M` = Mask
- `S` = Spectrum
- `T` = Threshold
- `W` = Workspace

**XX - Running Number (3rd-4th characters):**
- `01`, `02`, `03`, ... `99`
- Sequential numbering within each category

**Z - Variation Letter (5th character):**
- `A`, `B`, `C`, ...
- Distinguishes variations of the same basic functionality

#### Examples:
- `LI01A` = Load-Initialize-01-A (Initialize logging)
- `SM02A` = Selecting-Mask-02-A (Create mask type 2, variant A)
- `PA01B` = Processing-Averaging-01-B (Averaging method 1, variant B)
- `VT03A` = Visualize-Topo-03-A (Topography visualization, variant A)

---

## Block Anatomy

Every block follows this standardized structure:

```matlab
%% ABXXZ Category-Subcategory-XX-Z; Short description
% Edited by [Author Name], [Date]
%
% [Detailed description of what the block does]
% [Usage notes, dependencies, special considerations]

% Presets:
dataset = 'dataset1';           % specify the dataset to be used
variableIn = 'inputVar';        % specify the input variable
variableOut = 'outputVar';      % specify the output variable
param1 = value1;                % parameter description

% Optional variable inputs (set values to [] if not used)
optParam1 = [];                 % optional parameter description
optParam2 = [];                 % optional parameter description

%%%%%%%%%%%%%%%%%% DO NOT EDIT BELOW %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% LOG data in/out:
LOGcomment = sprintf("DataIn: %s.%s; dataOut: %s.%s", dataset, variableIn, dataset, variableOut);
LOGcomment = logUsedBlocks(LOGpath, LOGfile, "ABXXZ", LOGcomment, 0);

% Function execution
[data.(dataset).(variableOut), LOGcomment] = functionName(data.(dataset).(variableIn), param1);
LOGcomment = logUsedBlocks(LOGpath, LOGfile, "  ^  ", LOGcomment, 0);

% Clear preset variables
clearvars dataset variableIn variableOut param1 optParam1 optParam2
```

### Key Sections:

1. **Header Block** (`%%`)
   - Block identifier (ABXXZ)
   - Full descriptive name
   - Short description
   - Edit history
   - Detailed documentation

2. **Presets Section**
   - User-configurable parameters
   - Clear variable names with descriptions
   - Required vs. optional parameters clearly marked

3. **Do Not Edit Section**
   - Logging statements
   - Function execution
   - Output handling
   - Variable cleanup

4. **Logging Pattern**
   - First log: Block identifier + input/output summary
   - Subsequent logs: "  ^  " symbol + function details
   - Final log: Saved outputs (if applicable)

5. **Variable Cleanup**
   - `clearvars` removes all preset variables
   - Prevents variable name conflicts between blocks
   - Keeps workspace clean

---

## Section Organization

### Section 1: Documentation
- Block logging guidelines
- Naming conventions
- Usage instructions
- Best practices

### Section 2: Block List (Table of Contents)
- Organized by category
- Quick reference for all available blocks
- Includes retired blocks for reference

### Section 3: Initialization Blocks
- `LI01A`: Initialize logging system
- Sets up workspace
- Creates data structure

### Section 4: Data Loading Blocks
- `LD01A`, `LD02A`: Load data from files
- `LW01A`: Load saved workspace
- Handles different file formats

### Section 5: Workspace Management
- `SW01A`: Save workspace
- `LW01A`: Load workspace
- Includes log file management

### Section 6: Selecting Blocks
- Masking operations
- Logic operations (AND, OR, NOT)
- Region selection
- Threshold operations

### Section 7: Processing Blocks
- Data transformation
- Filtering and smoothing
- Derivative calculations
- Corrections and normalization

### Section 8: Visualizing Blocks
- 2D image plotting
- Spectrum/profile plotting
- Interactive viewers
- Video export

### Section 9: Custom Blocks
- Project-specific functionality
- Follows same structure as standard blocks

### Section 10: Retired Blocks
- Deprecated functionality
- Kept for reference
- Includes migration notes

---

## Data Structure Convention

### Primary Data Container: `data`

All data is stored in a structured format:

```matlab
data.<dataset>.<variable>
```

**Examples:**
```matlab
data.grid.I              % Current-voltage data
data.grid.V              % Voltage axis
data.grid.dIdV           % Derivative data
data.grid.mask           % Mask array
data.topo.z              % Topography height data
data.topo.x              % X-axis coordinates
```

### Benefits:
- Organized namespace
- Multiple datasets can coexist
- Clear data provenance
- Easy to save/load entire workspace

---

## Logging System

### Core Functions:

1. **`setLogFile()`**
   - Initializes logging system
   - Returns LOGpath and LOGfile
   - Creates log file

2. **`logUsedBlocks(LOGpath, LOGfile, blockID, comment, flag)`**
   - Records block execution
   - Timestamps all operations
   - Parameters:
     - `blockID`: Block identifier (e.g., "PA01A") or "  ^  " for continuation
     - `comment`: Description of operation
     - `flag`: 1 for new section, 0 for continuation

3. **`saveUsedBlocksLog(LOGpath, LOGfile, targetPath, filename)`**
   - Saves copy of log with output files
   - Links outputs to their creation log
   - Enables full reproducibility

### Logging Pattern:

```matlab
% First log entry for a block
LOGcomment = sprintf("DataIn: ...; DataOut: ...");
LOGcomment = logUsedBlocks(LOGpath, LOGfile, "PA01A", LOGcomment, 0);

% Subsequent entries use "  ^  "
LOGcomment = sprintf("Function: functionName(param1=%s, param2=%s)", val1, val2);
LOGcomment = logUsedBlocks(LOGpath, LOGfile, "  ^  ", LOGcomment, 0);

% Log saved outputs
LOGcomment = sprintf("Figure saved as: %s/%s.fig", path, name);
LOGcomment = logUsedBlocks(LOGpath, LOGfile, "  ^  ", LOGcomment, 0);
```

---

## Helper Functions

### Required Functions:

1. **`optionalStructCall(data, dataset, variable)`**
   - Safely accesses optional variables
   - Returns empty array if variable doesn't exist
   - Prevents errors from missing optional inputs

2. **`requiredStructCall(data, dataset, variable)`**
   - Accesses required variables
   - Throws error if variable doesn't exist
   - Ensures data integrity

3. **`uniqueNamePrompt(defaultName, suffix, path)`**
   - Prompts user for filename
   - Checks for existing files
   - Appends running number if needed
   - Prevents accidental overwrites

4. **`saveUsedBlocksLog(LOGpath, LOGfile, targetPath, filename)`**
   - Copies log file to output location
   - Creates `<filename>_LOGfile.txt`
   - Links outputs to their creation log

---

## Best Practices

### 1. Block Design
- **Single Responsibility**: Each block does one thing well
- **Clear Naming**: Block ID and description match functionality
- **Comprehensive Logging**: Log all inputs, outputs, and parameters
- **Clean Workspace**: Always clear preset variables at end

### 2. Variable Naming
- **Descriptive Names**: `variableIn1`, `variableOut1` over `x`, `y`
- **Consistent Patterns**: Same names for similar concepts across blocks
- **Dataset Specification**: Always specify which dataset is being used

### 3. Optional Parameters
- **Use [] for Unused**: Set optional parameters to `[]` when not needed
- **Check Before Use**: Use `isempty()` or `optionalStructCall()`
- **Document Clearly**: Explain what each optional parameter does

### 4. Plotting and Output
- **Always Offer Save Option**: Use `saveplots` boolean flag
- **Unique Names**: Use `uniqueNamePrompt()` for all saved outputs
- **Link Logs**: Call `saveUsedBlocksLog()` for all saved files
- **User Control**: Let user specify save location when appropriate

### 5. Error Handling
- **Validate Inputs**: Check data dimensions and types
- **Informative Messages**: Use `fprintf()` and `disp()` for user feedback
- **Graceful Degradation**: Handle missing optional inputs elegantly

### 6. Documentation
- **Block Headers**: Explain what, why, and how
- **Parameter Descriptions**: Comment every preset variable
- **Edit History**: Track who changed what and when
- **Usage Examples**: Include in comments when helpful

---

## Adapting the Template

### For a New Project:

1. **Copy the Template**
   - Start with `TEMPLATE_mainTrunk_structure.m`
   - Rename to match your project (e.g., `mainAnalysis.m`)

2. **Customize Section 1**
   - Update documentation for your specific needs
   - Add project-specific conventions
   - Define your category/subcategory codes

3. **Implement Required Functions**
   - Create your logging system (`setLogFile`, `logUsedBlocks`)
   - Implement helper functions (`optionalStructCall`, etc.)
   - Add project-specific utilities

4. **Create Your Blocks**
   - Follow the standard block structure
   - Use consistent naming conventions
   - Maintain comprehensive logging

5. **Update Block List**
   - Keep Section 2 current as you add blocks
   - Organize by category
   - Document retired blocks

### Adding a New Block:

1. **Choose Block ID**: Follow `ABXXZ` convention
2. **Copy Block Template**: Use existing block as starting point
3. **Update Header**: Change ID, name, description
4. **Define Presets**: Add your parameters
5. **Implement Function Call**: Call your processing function
6. **Add Logging**: Log inputs, outputs, parameters
7. **Clean Up**: Clear all preset variables
8. **Update Block List**: Add entry to Section 2

### Creating Block Variations:

When you need a variation of an existing block:
- Keep the same `ABXX` prefix
- Change the `Z` suffix (A → B → C)
- Document differences in header
- Example: `SM01A` → `SM01B` (different masking method)

---

## Example Workflow

### Typical Analysis Session:

```matlab
%% 1. Initialize
% Run LI01A to set up logging

%% 2. Load Data
% Run LD01A to load primary dataset
% Run LD02A to load secondary dataset (if needed)

%% 3. Process Data
% Run PA01A to smooth data
% Run PD01A to calculate derivative
% Run PF01A to flatten/normalize

%% 4. Select Regions
% Run SM01A to create mask
% Run SM02A to create second mask
% Run SL01A to combine masks (if needed)

%% 5. Visualize
% Run VT01A to plot 2D images
% Run VS01A to plot spectra
% Run VG01A for interactive viewing

%% 6. Save Work
% Run SW01A to save workspace
```

Each block execution is logged, creating a complete record of the analysis pipeline.

---

## Troubleshooting

### Common Issues:

**Issue**: Variable not found error
- **Solution**: Check that previous blocks have been run
- **Solution**: Verify dataset and variable names match

**Issue**: Log file not created
- **Solution**: Ensure LI01A has been run first
- **Solution**: Check write permissions for LOGpath

**Issue**: Figure not saving
- **Solution**: Verify saveplots = true
- **Solution**: Check target folder exists and is writable

**Issue**: Workspace load fails
- **Solution**: Ensure matching _LOGfile.txt exists
- **Solution**: Check that workspace contains data, LOGpath, LOGfile variables

---

## Advanced Features

### Conditional Execution:

```matlab
if isempty(optionalParam)
    % Execute without optional parameter
    result = function(data, requiredParam);
else
    % Execute with optional parameter
    result = function(data, requiredParam, optionalParam);
end
```

### Multiple Outputs:

```matlab
[data.(dataset).(varOut1), data.(dataset).(varOut2), LOGcomment] = ...
    function(data.(dataset).(varIn1), data.(dataset).(varIn2));
```

### Nested Data Access:

```matlab
% Using helper functions
xAxis = optionalStructCall(data, dataset, 'xAxis');
yData = requiredStructCall(data, dataset, 'yData');
```

### Interactive vs. Programmatic:

```matlab
% Interactive mode
if isempty(pointsList)
    pointsList = ginput(n);  % User clicks points
end

% Use pointsList (either from user or preset)
result = processPoints(data, pointsList);
```

---

## Version Control Recommendations

### What to Track:
- Main script file (e.g., `mainAnalysis.m`)
- Function files in `/functions/` directory
- Template and documentation files
- Example configuration files

### What NOT to Track:
- Log files (`*_LOGfile.txt`)
- Saved workspaces (`*.mat`)
- Generated figures (`*.fig`)
- Output data files

### Commit Messages:
- Use block IDs in commit messages
- Example: "Add SM03A: Polygon mask selection"
- Example: "Update PA01A: Add optional smoothing parameter"

---

## Conclusion

This template provides a robust, scalable framework for data analysis workflows. By following the standardized structure, you ensure:

- **Reproducibility**: Complete log of all operations
- **Maintainability**: Consistent structure across all blocks
- **Scalability**: Easy to add new functionality
- **Collaboration**: Clear conventions for team members
- **Documentation**: Self-documenting code structure

Adapt the template to your specific needs while maintaining the core principles of modular design, comprehensive logging, and standardized structure.

