function [data, params, isolation_results] = findIsolatedPoints(log, data, params, varargin)
%FINDISOLATEDPOINTS Find most isolated defect points for kernel initialization
%
%   Identifies the most isolated activation point for each kernel to use as
%   centers for 3D kernel initialization across all slices.
%
%   [data, params, isolation_results] = findIsolatedPoints(log, data, params, ...)
%
%   INPUTS:
%       log                 - Log struct with .path and .file fields
%       data                - Data struct from previous blocks (must contain data.slice from DS01A; or flat after extract)
%       params              - Parameter struct from previous blocks (hierarchical or flat)
%
%       OPTIONAL (Name-Value pairs):
%       'isolation_threshold_factor' - Threshold = max/factor for defect detection (default: 10)
%       'target_kernel_size_type' - Kernel sizing strategy (default: 'kernel_sizes_all')
%                                   Options:
%                                   'ref_kernel_sizes' - Use reference slice sizes for all
%                                   'kernel_sizes_cap' - Use maximum size across slices
%                                   'kernel_sizes_all' - Use slice-specific sizes
%       'show_distributions' - Show activation value distributions (default: true)
%       'show_isolation_maps' - Show isolation analysis visualizations (default: true)
%
%   OUTPUTS:
%       data                - Updated data struct with new fields in mcsbd_slice:
%                             data.slice.kernel_centers - [K×2] most isolated points
%                             data.slice.A0_used - Ground truth kernels (padded)
%                             data.slice.defect_positions - All defect locations
%                             data.slice.most_isolated_points - Selected centers
%                             data.slice.isolation_scores - Distance scores
%                             data.slice.num_defects - Count per kernel
%                             data.slice.offset - Alignment offset
%       params              - Updated parameter struct with:
%                             params.target_kernel_sizes - Sizes for 3D initialization
%
%   DESCRIPTION:
%       This wrapper finds the most spatially isolated activation point for
%       each kernel to serve as initialization centers for 3D processing:
%       1. Detects defects above threshold in reference activation maps
%       2. Calculates isolation scores (distance to other-kernel defects)
%       3. Filters boundary points and selects most isolated point per kernel
%       4. Aligns with ground truth and prepares kernel sizes for 3D
%
%   EXAMPLE:
%       % Standard isolation analysis
%       [data, params, isolation_results] = findIsolatedPoints(log, data, params);
%
%       % Custom threshold with slice-specific sizes
%       [data, params, isolation_results] = findIsolatedPoints(log, data, params, ...
%           'isolation_threshold_factor', 15, ...
%           'target_kernel_size_type', 'kernel_sizes_all');
%
%   See also: alignActivationMaps, padKernels

    % Parse input arguments
    p = inputParser;
    addRequired(p, 'log', @isstruct);
    addRequired(p, 'data', @isstruct);
    addRequired(p, 'params', @isstruct);
    addParameter(p, 'isolation_threshold_factor', 10, @isnumeric);
    addParameter(p, 'target_kernel_size_type', 'kernel_sizes_all', ...
        @(x) ismember(x, {'ref_kernel_sizes', 'kernel_sizes_cap', 'kernel_sizes_all'}));
    addParameter(p, 'show_distributions', true, @islogical);
    addParameter(p, 'show_isolation_maps', true, @islogical);
    parse(p, log, data, params, varargin{:});
    
    % Extract parameters from hierarchical structure (if needed)
    if isfield(params, 'synGen')
        params = organizeParams(params, 'extract');
    end
    
    % Extract data from hierarchical structure (if needed)
    if isfield(data, 'synGen') || isfield(data, 'slice') || isfield(data, 'mcsbd_slice')
        data = organizeData(data, 'extract');
    end
    
    % Extract parameters
    isolation_threshold_factor = p.Results.isolation_threshold_factor;
    target_kernel_size_type = p.Results.target_kernel_size_type;
    show_distributions = p.Results.show_distributions;
    show_isolation_maps = p.Results.show_isolation_maps;
    
    % Validate required fields in data struct (flat structure, after unpacking)
    if ~isfield(data, 'X')
        error('Data struct must contain field: X (reference slice activation maps from DS01A)');
    end
    if ~isfield(data, 'Y_ref')
        error('Data struct must contain field: Y_ref');
    end
    if ~isfield(data, 'X0_ref')
        error('Data struct must contain field: X0_ref');
    end
    
    % Validate required fields in params struct (flat structure, after unpacking)
    required_params = {'num_kernels', 'ref_slice', 'kernel_sizes'};
    for i = 1:length(required_params)
        if ~isfield(params, required_params{i})
            error('Params struct must contain field: %s', required_params{i});
        end
    end
    
    % LOG: Block start
    LOGcomment = sprintf("IS01A: threshold_factor=%g, sizing=%s", isolation_threshold_factor, target_kernel_size_type);
    LOGcomment = logUsedBlocks(log.path, log.file, "IS01A", LOGcomment, 0);
    
    fprintf('  Calculating isolation scores and finding most isolated points...\n');
    
    % Determine target kernel sizes for 3D initialization
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
    
    % Analyze activation value distributions
    if show_distributions
        figure('Name', 'IS01A: Activation Value Distributions');
    end
    
    for k = 1:params.num_kernels
        if show_distributions
            % Histogram of activation values
            subplot(2, params.num_kernels, k);
            activation_values = data.X(:,:,k);
            histogram(activation_values(activation_values > 0), 50);
            set(gca, 'YScale', 'log');
            title(sprintf('Kernel %d Distribution', k));
            xlabel('Activation Value');
            ylabel('Frequency (log)');
            
            % Threshold line
            threshold = max(data.X(:,:,k), [], 'all') / isolation_threshold_factor;
            hold on;
            xline(threshold, 'r--', 'Threshold');
            hold off;
            
            % Cumulative distribution
            subplot(2, params.num_kernels, k + params.num_kernels);
            [counts, edges] = histcounts(activation_values(activation_values > 0), 50, 'Normalization', 'cdf');
            stairs(edges(1:end-1), counts);
            title(sprintf('Kernel %d CDF', k));
            xlabel('Activation Value');
            ylabel('Cumulative Frequency');
            hold on;
            xline(threshold, 'r--', 'Threshold');
            hold off;
        else
            threshold = max(data.X(:,:,k), [], 'all') / isolation_threshold_factor;
        end
        
        % Get defect positions above threshold
        [rows, cols] = find(data.X(:,:,k) > threshold);
        defect_positions{k} = [rows, cols];
        num_defects(k) = size(defect_positions{k}, 1);
        fprintf('    Kernel %d: %d defects above threshold %.4f\n', k, num_defects(k), threshold);
    end
    
    % Calculate isolation scores
    for k = 1:params.num_kernels
        if num_defects(k) == 0
            warning('No defects found for kernel %d', k);
            continue;
        end
        
        % Sum activation maps of all other kernels
        X_others = zeros(size(data.X(:,:,1)));
        for l = 1:params.num_kernels
            if l ~= k
                X_others = X_others + data.X(:,:,l);
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
            
            if y <= half_kernel_size(1) || y >= size(data.X,1) - half_kernel_size(1) || ...
               x <= half_kernel_size(2) || x >= size(data.X,2) - half_kernel_size(2)
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
    
    % Visualize isolation analysis
    if show_isolation_maps
        figure('Name', 'IS01A: Isolation Analysis');
        for k = 1:params.num_kernels
            % Activation map with all defects and most isolated point
            subplot(2, params.num_kernels, k);
            imagesc(data.X(:,:,k));
            hold on;
            scatter(defect_positions{k}(:,2), defect_positions{k}(:,1), 50, 'w', 'o');
            if ~isempty(most_isolated_points{k})
                scatter(most_isolated_points{k}(2), most_isolated_points{k}(1), 100, 'r', '*', 'LineWidth', 2);
            end
            title(sprintf('Kernel %d', k));
            colorbar;
            axis square;
            colormap(gray);
            hold off;
            
            % Other kernels with most isolated point marked
            subplot(2, params.num_kernels, k + params.num_kernels);
            X_others = zeros(size(data.X(:,:,1)));
            for l = 1:params.num_kernels
                if l ~= k
                    X_others = X_others + data.X(:,:,l);
                end
            end
            imagesc(X_others);
            hold on;
            if ~isempty(most_isolated_points{k})
                scatter(most_isolated_points{k}(2), most_isolated_points{k}(1), 100, 'r', '*', 'LineWidth', 2);
            end
            title(sprintf('Other Kernels (not %d)', k));
            colorbar;
            axis square;
            colormap(gray);
            hold off;
        end
        
        % Display most isolated points on observation
        figure('Name', 'IS01A: Isolated Points on Observation');
        imagesc(data.Y_ref);
        colormap(gray);
        hold on;
        for k = 1:params.num_kernels
            scatter(most_isolated_points{k}(2), most_isolated_points{k}(1), 100, 'r', '*', 'LineWidth', 2);
            text(most_isolated_points{k}(2)+5, most_isolated_points{k}(1), sprintf('K%d', k), ...
                'Color', 'red', 'FontWeight', 'bold', 'FontSize', 12);
        end
        hold off;
        title('Most Isolated Points for Each Kernel');
        axis square;
        colorbar;
    end
    
    % Prepare ground truth kernels with appropriate sizes
    if strcmp(target_kernel_size_type, 'kernel_sizes_all')
        A0_used = data.A0;
    else
        A0_used = padKernels(data.A0_noiseless, params.SNR, target_kernel_sizes);
    end
    
    % Align most isolated points with ground truth
    kernel_sizes_ref_temp = reshape(params.kernel_sizes(params.ref_slice,:,:), [params.num_kernels, 2]);
    [~, offset, ~] = alignActivationMaps(data.X0_ref, data.X, kernel_sizes_ref_temp);
    used_most_isolated_points = cell(1, params.num_kernels);
    for k = 1:params.num_kernels
        used_most_isolated_points{k} = most_isolated_points{k} + offset(k,:);
    end
    
    % Convert to matrix format for kernel centers
    kernel_centers = zeros(params.num_kernels, 2);
    for k = 1:params.num_kernels
        if ~isempty(most_isolated_points{k})
            kernel_centers(k,:) = most_isolated_points{k};
        else
            error('No isolated point found for kernel %d', k);
        end
    end
    
    fprintf('  Isolation analysis complete.\n');
    
    % Store new results in flat structure (will be repacked to hierarchical)
    data.kernel_centers = kernel_centers;
    data.A0_used = A0_used;
    data.defect_positions = defect_positions;
    data.most_isolated_points = most_isolated_points;
    data.used_most_isolated_points = used_most_isolated_points;
    data.isolation_scores = isolation_scores;
    data.num_defects = num_defects;
    data.offset = offset;
    
    % Repack data to hierarchical structure for storage
    data = organizeData(data, 'write');
    
    % Store in params struct
    params.target_kernel_sizes = target_kernel_sizes;
    params.isolation_threshold_factor = isolation_threshold_factor;
    params.target_kernel_size_type = target_kernel_size_type;
    
    fprintf('  Found %d isolated points for 3D initialization.\n', params.num_kernels);
    
    % LOG: Isolation results
    LOGcomment = sprintf("Found %d isolated points, Defects per kernel: %s", ...
        params.num_kernels, mat2str(num_defects));
    LOGcomment = logUsedBlocks(log.path, log.file, "  ^  ", LOGcomment, 0);
    
    % Return isolation_results struct for backward compatibility
    isolation_results = struct();
    isolation_results.kernel_centers = kernel_centers;
    isolation_results.num_defects = num_defects;
    isolation_results.isolation_scores = isolation_scores;
    
    % Convert params to hierarchical structure for storage
    params = organizeParams(params, 'write');
end

