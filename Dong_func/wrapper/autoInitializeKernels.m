function [data, params] = autoInitializeKernels(log, data, params, varargin)
%AUTOINITIALIZEKERNELS Automatically initialize kernels from ground truth activations
%   Now includes internal logging
%
%   Automatically finds most isolated points in ground truth activation maps
%   and initializes kernels using those positions. Does NOT require running
%   MT-SBD algorithm first - works directly on synthetic data.
%
%   [data, params] = autoInitializeKernels(log, data, params, ...)
%
%   INPUTS:
%       log                 - Log struct for tracking execution
%       data                - Data struct with synGen fields
%       params              - Parameter struct from generateSyntheticData
%
%       OPTIONAL (Name-Value pairs):
%       'isolation_threshold_factor' - Threshold = max/factor for defect detection (default: 10)
%       'target_kernel_size_type' - Kernel sizing strategy (default: 'kernel_sizes_all')
%                                   Options:
%                                   'ref_kernel_sizes' - Use reference slice sizes for all
%                                   'kernel_sizes_cap' - Use maximum size across slices
%                                   'kernel_sizes_all' - Use slice-specific sizes
%       'window_type'       - Window function type (default: 'gaussian')
%       'window_sigma'      - Window parameter for Gaussian (default: 2.5)
%       'show_isolation'    - Show isolation analysis (default: false)
%
%   OUTPUTS:
%       data                - Updated data struct with new field:
%                             data.initialization.A_init - {1×K} initialized kernels
%                             data.initialization.kernel_centers - [K×2] positions  
%       params              - Updated parameter struct with:
%                             params.initialization.init_method - 'auto'
%                             params.initialization.init_window - {type, sigma}
%                             params.initialization.kernel_sizes - sizes for initialization
%
%   DESCRIPTION:
%       This wrapper automatically initializes kernels for the reference slice
%       by finding the most isolated activation points in ground truth data:
%       1. Detects defects above threshold in ground truth activation maps
%       2. Calculates isolation scores (distance to other-kernel defects)
%       3. Selects most isolated point per kernel
%       4. Initializes kernels at those positions with Gaussian window
%
%   EXAMPLE:
%       % Auto initialization with defaults
%       [data, params] = autoInitializeKernels(log, data, params);
%
%       % Custom window and threshold
%       [data, params] = autoInitializeKernels(log, data, params, ...
%           'window_sigma', 3.0, 'isolation_threshold_factor', 15);
%
%   See also: findIsolatedPoints, initialize_kernels, apply_window

    % Parse input arguments
    p = inputParser;
    addRequired(p, 'log', @isstruct);
    addRequired(p, 'data', @isstruct);
    addRequired(p, 'params', @isstruct);
    addParameter(p, 'isolation_threshold_factor', 10, @isnumeric);
    addParameter(p, 'target_kernel_size_type', 'kernel_sizes_all', ...
        @(x) ismember(x, {'ref_kernel_sizes', 'kernel_sizes_cap', 'kernel_sizes_all'}));
    addParameter(p, 'window_type', 'gaussian', @ischar);
    addParameter(p, 'window_sigma', 2.5, @isnumeric);
    addParameter(p, 'show_isolation', false, @islogical);
    parse(p, log, data, params, varargin{:});
    
    % Extract parameters from hierarchical structure (if needed)
    if isfield(params, 'initialization')
        params = organizeParams(params, 'extract');
    end
    
    % Extract data from hierarchical structure (if needed)
    if isfield(data, 'synGen') || isfield(data, 'initialization') || isfield(data, 'mcsbd_slice')
        data = organizeData(data, 'extract');
    end
    
    % Extract parameters
    isolation_threshold_factor = p.Results.isolation_threshold_factor;
    target_kernel_size_type = p.Results.target_kernel_size_type;
    window_type = p.Results.window_type;
    window_sigma = p.Results.window_sigma;
    show_isolation = p.Results.show_isolation;
    
    fprintf('  Auto-initializing kernels from ground truth activations...\n');
    
    % Validate required fields in data struct (flat structure)
    required_fields = {'X0_ref', 'A0_ref', 'Y_ref'};
    for i = 1:length(required_fields)
        if ~isfield(data, required_fields{i})
            error('Data struct must contain field: %s', required_fields{i});
        end
    end
    
    % Use ground truth activations directly
    X_gt = data.X0_ref;
    
    % Determine target kernel sizes for initialization
    switch target_kernel_size_type
        case 'ref_kernel_sizes'
            target_kernel_sizes = squeeze(params.kernel_sizes(params.ref_slice,:,:));
        case 'kernel_sizes_cap'
            target_kernel_sizes = squeeze(max(params.kernel_sizes,[],1));
        case 'kernel_sizes_all'
            target_kernel_sizes = params.kernel_sizes;
    end
    
    % Initialize storage
    most_isolated_points = cell(1, params.num_kernels);
    isolation_scores = cell(1, params.num_kernels);
    defect_positions = cell(1, params.num_kernels);
    num_defects = zeros(1, params.num_kernels);
    
    % Find isolated points for each kernel
    for k = 1:params.num_kernels
        threshold = max(X_gt(:,:,k), [], 'all') / isolation_threshold_factor;
        
        % Get defect positions above threshold
        [rows, cols] = find(X_gt(:,:,k) > threshold);
        defect_positions{k} = [rows, cols];
        num_defects(k) = size(defect_positions{k}, 1);
        fprintf('    Kernel %d: %d defects above threshold %.4f\n', k, num_defects(k), threshold);
        
        if num_defects(k) == 0
            warning('No defects found for kernel %d', k);
            continue;
        end
    end
    
    % Calculate isolation scores
    for k = 1:params.num_kernels
        if num_defects(k) == 0
            continue;
        end
        
        % Sum activation maps of all other kernels
        X_others = zeros(size(X_gt(:,:,1)));
        for l = 1:params.num_kernels
            if l ~= k
                X_others = X_others + X_gt(:,:,l);
            end
        end
        
        % Get positions in other kernels
        [other_rows, other_cols] = find(X_others > max(X_others,[],'all') / isolation_threshold_factor);
        other_positions = [other_rows, other_cols];
        
        % Boundary check based on target kernel size
        if strcmp(target_kernel_size_type, 'kernel_sizes_all')
            half_kernel_size = floor(squeeze(target_kernel_sizes(params.ref_slice,k,:))' / 2);
        else
            half_kernel_size = floor(target_kernel_sizes(k,:) / 2);
        end
        
        % Filter out points too close to boundaries
        valid_points = true(num_defects(k), 1);
        for i = 1:num_defects(k)
            y = defect_positions{k}(i,1);
            x = defect_positions{k}(i,2);
            
            if y <= half_kernel_size(1) || y >= size(X_gt,1) - half_kernel_size(1) || ...
               x <= half_kernel_size(2) || x >= size(X_gt,2) - half_kernel_size(2)
                valid_points(i) = false;
            end
        end
        
        % Process valid points
        valid_defects = defect_positions{k}(valid_points,:);
        if isempty(valid_defects)
            error('No valid isolated points for kernel %d - all too close to boundaries', k);
        end
        
        % Calculate isolation scores (distance to nearest other-kernel defect)
        S_k = zeros(size(valid_defects, 1), 1);
        for i = 1:size(valid_defects, 1)
            diffs = other_positions - valid_defects(i,:);
            distances = sum(diffs.^2, 2);
            S_k(i) = min(distances);
        end
        
        % Find most isolated point
        [max_score, max_idx] = max(S_k);
        valid_indices = find(valid_points);
        max_idx = valid_indices(max_idx);
        
        most_isolated_points{k} = defect_positions{k}(max_idx,:);
        isolation_scores{k} = S_k;
        
        fprintf('    Kernel %d: Most isolated at (%d,%d), score=%.2f\n', ...
            k, most_isolated_points{k}(1), most_isolated_points{k}(2), max_score);
    end
    
    % Visualize isolation analysis if requested
    if show_isolation
        figure('Name', 'Auto Kernel Initialization: Isolation Analysis');
        for k = 1:params.num_kernels
            subplot(2, params.num_kernels, k);
            imagesc(X_gt(:,:,k));
            hold on;
            if ~isempty(defect_positions{k})
                scatter(defect_positions{k}(:,2), defect_positions{k}(:,1), 50, 'w', 'o');
            end
            if ~isempty(most_isolated_points{k})
                scatter(most_isolated_points{k}(2), most_isolated_points{k}(1), 100, 'r', '*', 'LineWidth', 2);
            end
            title(sprintf('Kernel %d GT Activation', k));
            colorbar;
            axis square;
            colormap(gray);
            hold off;
            
            % Show on observation
            subplot(2, params.num_kernels, k + params.num_kernels);
            imagesc(data.Y_ref);
            hold on;
            if ~isempty(most_isolated_points{k})
                scatter(most_isolated_points{k}(2), most_isolated_points{k}(1), 100, 'r', '*', 'LineWidth', 2);
                text(most_isolated_points{k}(2)+5, most_isolated_points{k}(1), sprintf('K%d', k), ...
                    'Color', 'red', 'FontWeight', 'bold', 'FontSize', 12);
            end
            title(sprintf('Position on Y_{ref}', k));
            colorbar;
            axis square;
            colormap(gray);
            hold off;
        end
        sgtitle('Auto Kernel Initialization');
    end
    
    % Initialize kernels at isolated positions
    A_init = cell(1, params.num_kernels);
    kernel_centers = zeros(params.num_kernels, 2);
    
    for k = 1:params.num_kernels
        if isempty(most_isolated_points{k})
            error('No isolated point found for kernel %d', k);
        end
        
        kernel_centers(k,:) = most_isolated_points{k};
        
        % Get kernel size for this kernel
        if strcmp(target_kernel_size_type, 'kernel_sizes_all')
            ksize = squeeze(target_kernel_sizes(params.ref_slice,k,:))';
        else
            ksize = target_kernel_sizes(k,:);
        end
        
        % Extract patch from observation at isolated position
        y_center = most_isolated_points{k}(1);
        x_center = most_isolated_points{k}(2);
        half_size = floor(ksize / 2);
        
        y_range = (y_center - half_size(1)):(y_center + half_size(1));
        x_range = (x_center - half_size(2)):(x_center + half_size(2));
        
        % Handle size adjustment
        if length(y_range) > ksize(1)
            y_range = y_range(1:ksize(1));
        end
        if length(x_range) > ksize(2)
            x_range = x_range(1:ksize(2));
        end
        
        kernel_patch = data.Y_ref(y_range, x_range);
        
        % Apply window function
        if strcmp(window_type, 'gaussian')
            % For gaussian window with sigma parameter
            A_init{k} = windowToKernel(kernel_patch, window_type, window_sigma);
        elseif ~isempty(window_type) && ~strcmp(window_type, 'none')
            % For other window types without parameters
            A_init{k} = windowToKernel(kernel_patch, window_type);
        else
            % No window
            A_init{k} = kernel_patch;
        end
        
        fprintf('    Kernel %d: Initialized with size [%d×%d]\n', k, size(A_init{k},1), size(A_init{k},2));
    end
    
    % Store results in flat structure (will convert to hierarchical at end)
    % Data: A_init at data.A_init (will be organized to data.initialization.A_init)
    data.A_init = A_init;
    % Store kernel centers in data (will be organized to data.initialization.kernel_centers)
    data.init_kernel_centers = kernel_centers;
    
    % Store parameters in flat structure (will convert to hierarchical at end)
    params.init_method = 'auto';
    params.init_window = {window_type, window_sigma};
    params.kernel_sizes = target_kernel_sizes;  % Note: initialization sizes (may differ from generation sizes)
    % Note: kernel_centers stored only in data.initialization, not in params
    
    % LOG: Auto initialization results
    LOGcomment = sprintf("Auto-initialized kernels: method=%s, num_kernels=%d" ...
        ,params.init_method, params.num_kernels);
    LOGcomment = logUsedBlocks(log.path, log.file, "  ^  ", LOGcomment, 0);
    
    % Convert data and params to hierarchical structure for storage
    data = organizeData(data, 'write');
    params = organizeParams(params, 'write');
    
    fprintf('  Auto kernel initialization complete.\n');
end

