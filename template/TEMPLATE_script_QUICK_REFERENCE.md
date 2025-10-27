# Quick Reference Guide: Block-Based Analysis Template

## Block ID Format: `ABXXZ`

| Position | Meaning | Examples |
|----------|---------|----------|
| **A** (Category) | L=Load, P=Process, V=Visualize, S=Select, W=Write | L, P, V, S, W |
| **B** (Subcategory) | A=Average, D=Data, F=Flatten, G=Grid, I=Init, M=Mask, S=Spectrum, T=Threshold, W=Workspace | A, D, F, G, I, M, S, T, W |
| **XX** (Number) | Sequential number within category | 01, 02, ..., 99 |
| **Z** (Variant) | Variation letter | A, B, C, ... |

### Example IDs:
- `LI01A` = Load-Initialize-01-A
- `SM02B` = Select-Mask-02-B  
- `PA01A` = Process-Average-01-A
- `VT03A` = Visualize-Topo-03-A

---

## Standard Block Structure

```matlab
%% ABXXZ Category-Subcategory-XX-Z; Description
% Edited by [Name], [Date]
% [Documentation]

% Presets:
dataset = 'name';
variableIn = 'input';
variableOut = 'output';
param = value;

% Optional (set to [] if unused):
optParam = [];

%%%%%%%%%%%%%%%%%% DO NOT EDIT BELOW %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% LOG
LOGcomment = sprintf("DataIn: %s.%s; DataOut: %s.%s", ...);
LOGcomment = logUsedBlocks(LOGpath, LOGfile, "ABXXZ", LOGcomment, 0);

% Execute
[data.(dataset).(variableOut), LOGcomment] = func(...);
LOGcomment = logUsedBlocks(LOGpath, LOGfile, "  ^  ", LOGcomment, 0);

% Cleanup
clearvars dataset variableIn variableOut param optParam
```

---

## Essential Functions

### Logging Functions
```matlab
% Initialize logging
[LOGpath, LOGfile] = setLogFile();

% Log block execution
LOGcomment = logUsedBlocks(LOGpath, LOGfile, "ABXXZ", comment, 0);

% Log continuation (use "  ^  " instead of block ID)
LOGcomment = logUsedBlocks(LOGpath, LOGfile, "  ^  ", comment, 0);

% Save log copy with output
saveUsedBlocksLog(LOGpath, LOGfile, targetPath, filename);
```

### Data Access Functions
```matlab
% Access optional variable (returns [] if missing)
var = optionalStructCall(data, dataset, 'varName');

% Access required variable (errors if missing)
var = requiredStructCall(data, dataset, 'varName');
```

### File Management
```matlab
% Get unique filename (auto-appends number if exists)
filename = uniqueNamePrompt("DefaultName", "", path);

% Save figure
savefig(strcat(path, "/", filename, ".fig"));
```

---

## Data Structure Convention

```matlab
data.<dataset>.<variable>
```

### Examples:
```matlab
data.grid.I              % Raw current data
data.grid.V              % Voltage axis
data.grid.dIdV           % Derivative
data.grid.mask           % Mask array
data.topo.z              % Height data
data.topo.x              % X coordinates
```

---

## Common Block Patterns

### Pattern 1: Simple Processing
```matlab
%% PA01A Processing-Average-01-A; Average data

% Presets:
dataset = 'grid';
variableIn = 'I';
variableOut = 'I_avg';
windowSize = 5;

%%%%%%%%%%%%%%%%%% DO NOT EDIT BELOW %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
LOGcomment = sprintf("DataIn: %s.%s; DataOut: %s.%s", dataset, variableIn, dataset, variableOut);
LOGcomment = logUsedBlocks(LOGpath, LOGfile, "PA01A", LOGcomment, 0);

[data.(dataset).(variableOut), LOGcomment] = smoothData(data.(dataset).(variableIn), windowSize);
LOGcomment = logUsedBlocks(LOGpath, LOGfile, "  ^  ", LOGcomment, 0);

clearvars dataset variableIn variableOut windowSize
```

### Pattern 2: Processing with Optional Parameters
```matlab
%% PF01A Processing-Flatten-01-A; Flatten data

% Presets:
dataset = 'topo';
variableIn1 = 'z';
variableIn2 = '';           % optional mask ('' if not used)
variableOut = 'z_flat';
nPoints = 100;

%%%%%%%%%%%%%%%%%% DO NOT EDIT BELOW %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
LOGcomment = sprintf("DataIn: %s.%s, %s; DataOut: %s.%s", dataset, variableIn1, variableIn2, dataset, variableOut);
LOGcomment = logUsedBlocks(LOGpath, LOGfile, "PF01A", LOGcomment, 0);

if isempty(variableIn2)
    [data.(dataset).(variableOut), LOGcomment] = flattenFunc(data.(dataset).(variableIn1), '', nPoints);
else
    [data.(dataset).(variableOut), LOGcomment] = flattenFunc(data.(dataset).(variableIn1), data.(dataset).(variableIn2), nPoints);
end
LOGcomment = logUsedBlocks(LOGpath, LOGfile, "  ^  ", LOGcomment, 0);

clearvars dataset variableIn1 variableIn2 variableOut nPoints
```

### Pattern 3: Visualization with Save Option
```matlab
%% VT01A Visualize-Topo-01-A; Plot 2D image

% Presets:
dataset = 'topo';
variableIn = 'z';
saveplots = true;           % true to save, false to skip

%%%%%%%%%%%%%%%%%% DO NOT EDIT BELOW %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
LOGcomment = sprintf("DataIn: %s.%s", dataset, variableIn);
LOGcomment = logUsedBlocks(LOGpath, LOGfile, "VT01A", LOGcomment, 0);

[fig, LOGcomment] = plotImage(data.(dataset).(variableIn));
LOGcomment = logUsedBlocks(LOGpath, LOGfile, "  ^  ", LOGcomment, 0);

if saveplots == true
    plot_name = uniqueNamePrompt("TopoImage", "", LOGpath);
    savefig(fig, strcat(LOGpath, "/", plot_name, ".fig"));
    saveUsedBlocksLog(LOGpath, LOGfile, LOGpath, plot_name);
    LOGcomment = sprintf("Saved as: %s/%s.fig", LOGpath, plot_name);
    LOGcomment = logUsedBlocks(LOGpath, LOGfile, "  ^  ", LOGcomment, 0);
end

clearvars dataset variableIn saveplots plot_name fig
```

### Pattern 4: Multiple Outputs
```matlab
%% PD01A Processing-Derivative-01-A; Calculate derivative

% Presets:
dataset = 'grid';
variableIn1 = 'I';
variableIn2 = 'V';
variableOut1 = 'dIdV';
variableOut2 = 'V_reduced';

%%%%%%%%%%%%%%%%%% DO NOT EDIT BELOW %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
LOGcomment = sprintf("DataIn: %s.%s, %s.%s; DataOut: %s.%s, %s.%s", ...
    dataset, variableIn1, dataset, variableIn2, dataset, variableOut1, dataset, variableOut2);
LOGcomment = logUsedBlocks(LOGpath, LOGfile, "PD01A", LOGcomment, 0);

[data.(dataset).(variableOut1), data.(dataset).(variableOut2), LOGcomment] = ...
    Derivative(data.(dataset).(variableIn1), data.(dataset).(variableIn2));
LOGcomment = logUsedBlocks(LOGpath, LOGfile, "  ^  ", LOGcomment, 0);

clearvars dataset variableIn1 variableIn2 variableOut1 variableOut2
```

---

## Typical Workflow

```matlab
%% 1. Initialize
% LI01A - Initialize logging

%% 2. Load Data  
% LD01A - Load primary data
% LD02A - Load secondary data (optional)

%% 3. Process
% PA01A - Smooth/average
% PD01A - Calculate derivative
% PF01A - Flatten/normalize

%% 4. Select Regions
% SM01A - Create mask
% SL01A - Combine masks (optional)

%% 5. Visualize
% VT01A - Plot images
% VS01A - Plot spectra

%% 6. Save
% SW01A - Save workspace
```

---

## Checklist for Creating a New Block

- [ ] Choose appropriate block ID (`ABXXZ`)
- [ ] Write clear header with description
- [ ] Define all preset variables with comments
- [ ] Mark optional parameters (set to `[]`)
- [ ] Add "DO NOT EDIT BELOW" separator
- [ ] Log inputs and outputs
- [ ] Execute function
- [ ] Log function execution with "  ^  "
- [ ] Handle optional save/plot functionality
- [ ] Clear all preset variables with `clearvars`
- [ ] Add block to Section 2 (Block List)
- [ ] Test with sample data

---

## Common Mistakes to Avoid

❌ **Don't**: Forget to clear preset variables
```matlab
% Missing clearvars - variables persist!
```

✅ **Do**: Always clear at end of block
```matlab
clearvars dataset variableIn variableOut param
```

---

❌ **Don't**: Use generic variable names in workspace
```matlab
x = data.grid.I;  % 'x' will conflict across blocks
```

✅ **Do**: Use structured data storage
```matlab
data.grid.I_smoothed = smoothData(data.grid.I);
```

---

❌ **Don't**: Skip logging
```matlab
result = processData(input);  % No log = not reproducible
```

✅ **Do**: Log everything
```matlab
[result, LOGcomment] = processData(input);
LOGcomment = logUsedBlocks(LOGpath, LOGfile, "  ^  ", LOGcomment, 0);
```

---

❌ **Don't**: Hardcode paths or filenames
```matlab
save('C:/Users/Me/data.mat');  % Won't work on other machines
```

✅ **Do**: Use relative paths and user prompts
```matlab
filename = uniqueNamePrompt("data", "", LOGpath);
save(strcat(LOGpath, '/', filename, '.mat'));
```

---

## Quick Syntax Reference

### String Formatting
```matlab
% sprintf for formatted strings
LOGcomment = sprintf("Value = %s, Number = %d", str, num);

% num2str for numbers to strings
LOGcomment = sprintf("Param = %s", num2str(value));

% mat2str for arrays to strings
LOGcomment = sprintf("Array = %s", mat2str([1 2 3]));
```

### Conditional Execution
```matlab
% Check if empty
if isempty(optParam)
    % Handle missing optional parameter
end

% Check if variable exists
if exist('varName', 'var')
    % Variable exists in workspace
end

% Check if field exists in struct
if isfield(data.grid, 'mask')
    % Field exists
end
```

### Array Operations
```matlab
% Get size
[nx, ny, nz] = size(data.grid.I);

% Get dimensions
dims = ndims(data.grid.I);  % Number of dimensions

% Reshape
data.grid.I_flat = reshape(data.grid.I, [], 1);  % Flatten to column
```

---

## File Organization

```
project/
├── mainAnalysis.m              # Main script (based on template)
├── functions/                  # Function library
│   ├── basic/                  # Basic utilities
│   ├── loading/                # Data loading functions
│   ├── processing/             # Processing functions
│   ├── selecting/              # Masking functions
│   └── visualization/          # Plotting functions
├── data/                       # Input data (not in git)
├── output/                     # Results (not in git)
├── logs/                       # Log files (not in git)
└── README.md                   # Project documentation
```

---

## Getting Help

### Template Files:
- `TEMPLATE_mainTrunk_structure.m` - Full template with examples
- `TEMPLATE_DOCUMENTATION.md` - Comprehensive guide
- `TEMPLATE_QUICK_REFERENCE.md` - This file

### Key Concepts:
- **Block**: Self-contained functional unit
- **Block ID**: Unique 5-character identifier (ABXXZ)
- **Preset**: User-configurable parameter
- **LOGcomment**: String describing operation for log
- **Dataset**: Named collection of related variables in data struct

### Common Questions:

**Q: How do I add a new block?**  
A: Copy an existing similar block, change the ID, update presets and function call.

**Q: What if my function doesn't return LOGcomment?**  
A: Create the LOGcomment string manually describing the operation.

**Q: Can I have multiple datasets?**  
A: Yes! Use different dataset names: `data.grid1`, `data.grid2`, `data.topo`, etc.

**Q: How do I handle errors?**  
A: Use try-catch blocks and log errors with `logUsedBlocks()`.

**Q: Can blocks call other blocks?**  
A: No - blocks should be independent. Share code via functions in `/functions/`.

---

## Version History

- v1.0 - Initial template based on UBC_LAIR mainTrunk.m
- Template created: 2025

---

**Remember**: Consistency is key. Follow the template structure, and your analysis will be reproducible, maintainable, and easy to understand!

