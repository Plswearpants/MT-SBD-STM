function [data, params] = initializeKernelsRef(log, data, params)
%INITIALIZEKERNELSREF Wrapper for initializing kernels from reference slice
%   Now includes internal logging
%
%   [data, params] = initializeKernelsRef(log, data, params)
%
%   INPUTS:
%       log                 - Log struct for tracking execution
%       data                - Data struct from previous block (hierarchical or flat)
%       params              - Parameter struct from previous block (hierarchical or flat)
%
%       REQUIRED (in params.synGen or params):
%       .num_kernels        - Number of kernels to initialize
%       .ref_slice          - Reference slice index
%       .kernel_sizes       - [num_slices × num_kernels × 2] kernel sizes
%
%       OPTIONAL (in params.synGen or params):
%       .kernel_selection_type - 'selected' (interactive) or 'random' (default: 'selected')
%       .window_type        - Window function to apply to kernels (default: {'gaussian', 2.5})
%                             Options: 'hann', 'hamming', 'blackman'
%                                      {'gaussian', alpha}
%                                      {'kaiser', beta}
%                                      {} or '' for no window
%
%   OUTPUTS:
%       data                - Updated data struct (hierarchical):
%                             data.initialization.A_init - {1×K} initialized kernels
%                             data.initialization.kernel_centers - [K×2] positions (empty for manual)
%       params              - Updated parameter struct (hierarchical):
%                             params.initialization.init_method - 'manualXX' with iteration number
%                             params.initialization.init_window - Window function {type, sigma}
%                             Note: kernel_centers stored only in data.initialization, not in params
%
%   DESCRIPTION:
%       This wrapper encapsulates the kernel initialization process for the
%       reference slice. It handles:
%       - Extracting kernel sizes for the reference slice
%       - Displaying ground truth kernels (only for 'selected' mode)
%       - Calling initialize_kernels() for interactive/random selection
%       - Displaying initialized kernels (always)
%       - Returning kernels with proper windowing applied
%
%   EXAMPLE:
%       % Set parameters in hierarchical structure
%       params.synGen.kernel_selection_type = 'selected';
%       params.synGen.window_type = {'gaussian', 2.5};
%       [data, params] = initializeKernelsRef(log, data, params);
%
%   See also: initialize_kernels, apply_window, proj2oblique, organizeData, organizeParams

    % Validate required inputs
    if ~isstruct(log) || ~isfield(log, 'path') || ~isfield(log, 'file')
        error('log must be a struct with .path and .file fields');
    end
    if ~isstruct(data)
        error('data must be a struct');
    end
    if ~isstruct(params)
        error('params must be a struct');
    end
    
    % Extract params from hierarchical structure (if needed)
    if isfield(params, 'synGen')
        params = organizeParams(params, 'extract');
    end
    
    % Extract data from hierarchical structure (if needed)
    if isfield(data, 'synGen') || isfield(data, 'mcsbd_slice')
        data = organizeData(data, 'extract');
    end
    
    % Read optional parameters with defaults
    if isfield(params, 'kernel_selection_type')
        kernel_selection_type = params.kernel_selection_type;
    else
        kernel_selection_type = 'selected';
    end
    
    if isfield(params, 'window_type')
        window_type = params.window_type;
    else
        window_type = {'gaussian', 2.5};
    end
    
    % Validate kernel_selection_type
    if ~ismember(kernel_selection_type, {'selected', 'random'})
        error('kernel_selection_type must be ''selected'' or ''random''');
    end
    
    % Automatic display logic:
    % - Show ground truth ONLY for 'selected' mode (to guide user selection)
    % - Always show initialized kernels (to verify initialization)
    show_ground_truth = strcmp(kernel_selection_type, 'selected');
    show_initialized = true;
    
    % Validate required fields in data struct (flat structure)
    required_fields = {'Y_ref', 'A0_ref'};
    for i = 1:length(required_fields)
        if ~isfield(data, required_fields{i})
            error('Data struct must contain field: %s', required_fields{i});
        end
    end
    
    % Validate required fields in params struct
    required_params = {'num_kernels', 'ref_slice', 'kernel_sizes'};
    for i = 1:length(required_params)
        if ~isfield(params, required_params{i})
            error('Params struct must contain field: %s', required_params{i});
        end
    end
    
    % Extract kernel sizes for reference slice
    % params.kernel_sizes is [num_slices × num_kernels × 2]
    kernel_sizes_ref = reshape(params.kernel_sizes(params.ref_slice,:,:), [params.num_kernels, 2]);
    
    % Initialize kernels
    fprintf('  Initializing kernels for reference slice...\n');
    if strcmp(kernel_selection_type, 'selected')
        fprintf('\n');
        fprintf('========================================================================\n');
        fprintf('  INTERACTIVE KERNEL SELECTION (Guided Mode)\n');
        fprintf('========================================================================\n');
        fprintf('  Layout:\n');
        fprintf('    LEFT:  Observation (Y_ref) - where you select the region\n');
        fprintf('    RIGHT: ALL Ground Truth kernels in a column\n');
        fprintf('           → Current kernel to select is HIGHLIGHTED\n');
        fprintf('\n');
        fprintf('  Instructions:\n');
        fprintf('    1. Look at the HIGHLIGHTED kernel on the right (yellow background)\n');
        fprintf('    2. Position the rectangle on the left to match that pattern\n');
        fprintf('    3. Double-click the rectangle to confirm and move to next kernel\n');
        fprintf('    4. Repeat until all kernels are selected\n');
        fprintf('========================================================================\n\n');
        
        % Use guided initialization with side-by-side GT display
        GT_kernels_ref = cell(1, params.num_kernels);
        for k = 1:params.num_kernels
            GT_kernels_ref{k} = data.A0{params.ref_slice, k};
        end
        
        A1 = initialize_kernels_guided(data.Y_ref, GT_kernels_ref, ...
                                       params.num_kernels, kernel_sizes_ref, window_type);
    else
        fprintf('  Generating %d random initialized kernels.\n', params.num_kernels);
        A1 = initialize_kernels(data.Y_ref, params.num_kernels, kernel_sizes_ref, ...
                                'random', window_type);
    end
    
    % Display initialized kernels
    if show_initialized
        fprintf('  Displaying initialized kernels...\n');
        figure('Name', 'IN01A: Initialized Kernels (Reference Slice)');
        for n = 1:params.num_kernels
            subplot(1, params.num_kernels, n);
            imagesc(A1{n});
            title(sprintf('Initial Kernel %d', n));
            colorbar;
            axis square;
            colormap(gray);
        end
        sgtitle(sprintf('Initialized Kernels (Reference Slice %d)', params.ref_slice));
    end
    
    fprintf('  Kernel initialization complete.\n');
    
    % Determine manual iteration number
    manual_iteration = 1;
    if isfield(params, 'init_method') && startsWith(params.init_method, 'manual')
        % Extract number from 'manual01', 'manual02', etc.
        prev_method = params.init_method;
        if length(prev_method) > 6
            manual_iteration = str2double(prev_method(7:end)) + 1;
        end
    end
    
    % Store results in flat structure (will convert to hierarchical at end)
    % Store A_init as primary location
    data.A_init = A1;
    % Store kernel centers (empty for manual, no automatic position detection)
    data.init_kernel_centers = [];
    
    % Store parameters according to documented structure (flat structure)
    params.init_method = sprintf('manual%02d', manual_iteration);
    params.init_window = window_type;
    % Note: kernel_centers stored only in data.initialization, not in params
    
    % Store additional parameters (not in docs but needed internally)
    params.kernel_selection_type = kernel_selection_type;
    params.kernel_sizes_ref = kernel_sizes_ref;
    
    % LOG: Manual initialization results
    LOGcomment = sprintf("Manual init complete: %s, kernels=%d", ...
        params.init_method, params.num_kernels);
    LOGcomment = logUsedBlocks(log.path, log.file, "  ^  ", LOGcomment, 0);
    
    % Convert data and params to hierarchical structure for storage
    data = organizeData(data, 'write');
    params = organizeParams(params, 'write');
    
    fprintf('  Manual kernel initialization stored as: %s\n', params.initialization.init_method);
end

