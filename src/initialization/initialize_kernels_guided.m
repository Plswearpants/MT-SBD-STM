function A1 = initialize_kernels_guided(Y, GT_kernels, num_kernels, kernel_sizes, window_type)
    %INITIALIZE_KERNELS_GUIDED Interactive kernel initialization with side-by-side GT reference
    %
    %   Improved UI for manual kernel selection that displays the ground truth
    %   kernel on the right side during selection to guide the user.
    %
    %   A1 = initialize_kernels_guided(Y, GT_kernels, num_kernels, kernel_sizes, window_type)
    %
    %   INPUTS:
    %       Y               - Reference slice observation image [H × W]
    %       GT_kernels      - Cell array of ground truth kernels {1 × num_kernels}
    %       num_kernels     - Number of kernels to initialize
    %       kernel_sizes    - Size of each kernel [num_kernels × 2] as [height, width]
    %       window_type     - (optional) Window function to apply
    %                         '' or [] for no window
    %                         'hann', 'hamming', 'blackman'
    %                         {'gaussian', alpha} or {'kaiser', beta}
    %
    %   OUTPUTS:
    %       A1              - Cell array of initialized kernels {1 × num_kernels}
    %
    %   DESCRIPTION:
    %       Displays Y (observation) on the left and the corresponding GT kernel
    %       on the right during selection. User selects regions sequentially,
    %       with visual guidance showing which kernel to select.
    %
    %   See also: initialize_kernels, initializeKernelsRef
    
    % Set default window_type to empty if not provided
    if nargin < 5
        window_type = '';
    end
    
    A1 = cell(1, num_kernels);
    [img_height, img_width] = size(Y);
    
    % Create figure with observation on left and all GT kernels on right
    fig = figure('Name', 'Manual Kernel Initialization', ...
                 'Position', [50, 50, 1600, 800]);
    
    % Store subplot handles for GT kernels (will add red box to highlight)
    gt_axes = cell(1, num_kernels);
    red_boxes = cell(1, num_kernels);
    
    % Set up the layout: LEFT = observation (takes 60% width), RIGHT = GT kernels column (40%)
    % Use custom positioning for better control
    
    for n = 1:num_kernels
        % Clear previous highlights and update layout
        
        % LEFT: Observation for selection (fixed position)
        ax_obs = subplot('Position', [0.05, 0.1, 0.5, 0.8]);
        imagesc(Y);
        colorbar;
        colormap(gray);
        axis square;
        title(sprintf('SELECT Kernel %d of %d (Position rectangle, double-click to confirm)', ...
              n, num_kernels), 'FontWeight', 'bold', 'FontSize', 14);
        xlabel('Position the rectangle over a region matching the HIGHLIGHTED kernel on right →', ...
               'FontSize', 11, 'FontWeight', 'bold');
        hold on;
        
        % RIGHT: Display ALL GT kernels in a column
        for k = 1:num_kernels
            % Calculate vertical position for each kernel
            % Leave space at top/bottom, distribute evenly
            kernel_height = 0.7 / num_kernels;  % Total height divided by number of kernels
            y_pos = 0.85 - (k * kernel_height);  % Start from top
            
            gt_axes{k} = subplot('Position', [0.65, y_pos, 0.28, kernel_height * 0.85]);
            imagesc(GT_kernels{k});
            colorbar;
            colormap(gray);
            axis square;
            
            % Title and appearance based on whether this is the current kernel
            if k == n
                % Current kernel to select - HIGHLIGHT
                title(sprintf('→ Kernel %d ← SELECT THIS ONE', k), ...
                      'FontWeight', 'bold', 'FontSize', 12, 'Color', [0.8, 0, 0], ...
                      'BackgroundColor', [1, 1, 0.8]);
            else
                % Other kernels - gray out
                title(sprintf('Kernel %d', k), 'FontSize', 10, 'Color', [0.5, 0.5, 0.5]);
                alpha(0.5);  % Make non-current kernels slightly transparent
            end
        end
        
        % Add overall figure title
        sgtitle(sprintf('Selecting Kernel %d of %d - Match the HIGHLIGHTED kernel on the right', ...
                n, num_kernels), 'FontSize', 16, 'FontWeight', 'bold', 'Color', [0.8, 0, 0]);
        
        % Go back to observation subplot for selection
        subplot(ax_obs);
        
        % Initial position (try to place in visible area)
        init_x = 10;
        init_y = 10;
        rect_width = min(kernel_sizes(n,2), img_width - 20);
        rect_height = min(kernel_sizes(n,1), img_height - 20);
        
        % Create rectangle
        h_rect = imrect(gca, [init_x init_y rect_width rect_height]);
        
        % Calculate center point of rectangle
        center_x = init_x + rect_width/2;
        center_y = init_y + rect_height/2;
        
        % Create center dot indicator
        h_dot = plot(center_x, center_y, 'r.', 'MarkerSize', 20);
        
        % Add listener to update center dot when rectangle moves
        addNewPositionCallback(h_rect, @(pos) updateCenterDot(pos, h_dot));
        
        % Set position constraints
        setPositionConstraintFcn(h_rect, @(pos) constrainPosition(pos, img_width, img_height, kernel_sizes(n,:)));
        
        % Display instruction
        fprintf('  Selecting Kernel %d/%d...\n', n, num_kernels);
        fprintf('    - All GT kernels shown on RIGHT (current one has yellow background)\n');
        fprintf('    - Position rectangle on LEFT to match the highlighted GT pattern\n');
        fprintf('    - Double-click rectangle to confirm selection\n');
        
        % Wait for user to finish positioning
        position = wait(h_rect);
        
        % Extract the selected region
        x1 = round(position(1));
        y1 = round(position(2));
        x2 = min(x1 + kernel_sizes(n,2) - 1, img_width);
        y2 = min(y1 + kernel_sizes(n,1) - 1, img_height);
        
        selected_kernel = Y(y1:y2, x1:x2);
        
        % Project onto the oblique manifold
        A1{n} = proj2oblique(selected_kernel);
        
        % Apply window if specified
        if ~isempty(window_type)
            if iscell(window_type)
                % For windows with parameters (gaussian, kaiser)
                A1{n} = windowToKernel(A1{n}, window_type{1}, window_type{2});
            else
                % For windows without parameters
                A1{n} = windowToKernel(A1{n}, window_type);
            end
        end
        
        fprintf('    ✓ Kernel %d selected and initialized\n\n', n);
    end
    
    % Close the selection figure
    close(fig);
    
    fprintf('  All %d kernels initialized successfully!\n', num_kernels);
end

%% Helper Functions

function updateCenterDot(pos, h_dot)
    % Update center dot position when rectangle moves
    center_x = pos(1) + pos(3)/2;
    center_y = pos(2) + pos(4)/2;
    set(h_dot, 'XData', center_x, 'YData', center_y);
end

function new_pos = constrainPosition(pos, img_width, img_height, kernel_size)
    % Constrain rectangle to stay within image bounds
    new_pos = pos;
    
    % Ensure rectangle doesn't go outside image
    if pos(1) < 1
        new_pos(1) = 1;
    end
    if pos(2) < 1
        new_pos(2) = 1;
    end
    if pos(1) + pos(3) > img_width
        new_pos(1) = img_width - pos(3);
    end
    if pos(2) + pos(4) > img_height
        new_pos(2) = img_height - pos(4);
    end
    
    % Maintain kernel size
    new_pos(3) = kernel_size(2);
    new_pos(4) = kernel_size(1);
end

