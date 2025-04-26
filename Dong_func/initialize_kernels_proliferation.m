function A1 = initialize_kernels_proliferation(Y, num_kernels, kernel_centers, window_type, target_kernel_size, varargin)
    % Initialize kernels for proliferation case using pre-defined centers
    % Inputs:
    %   Y: input image
    %   num_kernels: number of kernels to initialize
    %   kernel_centers: [num_kernels x 2] matrix of pre-defined centers [y,x]
    %   window_type: type of window to apply
    %               'hann', 'hamming', 'blackman', 'gaussian', 'kaiser'
    %               For gaussian/kaiser, use cell array: {'gaussian', alpha} or {'kaiser', beta}
    %   target_kernel_size: [num_kernels x 2] matrix of target kernel sizes [height, width]
    %   varargin: optional parameters
    %       'interactive': boolean, whether to allow interactive kernel size adjustment
    %                      default: true
    % Outputs:
    %   A1: cell array of initialized kernels
    %       Each cell contains a kernel of size [h w] projected onto oblique manifold
    %       Window function is applied to each kernel
    
    % Parse optional parameters
    p = inputParser;
    addParameter(p, 'interactive', true, @islogical);
    parse(p, varargin{:});
    interactive = p.Results.interactive;
    
    % Validate inputs
    if nargin < 5
        error('First 5 inputs are required');
    end
    
    if isempty(Y)
        error('Input image Y cannot be empty');
    end
    
    if ~ismatrix(Y)
        error('Input image Y must be a 2D matrix');
    end
    
    if num_kernels <= 0
        error('Number of kernels must be positive');
    end
    
    if size(kernel_centers,1) ~= num_kernels || size(kernel_centers,2) ~= 2
        error('kernel_centers must be a [num_kernels x 2] matrix');
    end
    
    if isempty(window_type)
        error('window_type is required');
    end
    
    if size(target_kernel_size,1) ~= num_kernels || size(target_kernel_size,2) ~= 2
        error('target_kernel_size must be a [num_kernels x 2] matrix');
    end
    
    % Initialize output
    A1 = cell(1, num_kernels);
    [img_height, img_width] = size(Y);
    
    % Store final kernel sizes
    kernel_sizes = zeros(num_kernels, 2);
    
    fprintf('Initializing kernels using pre-defined centers...\n');
    for n = 1:num_kernels
        % Get center coordinates
        center_y = kernel_centers(n,1);
        center_x = kernel_centers(n,2);
        
        % Validate center coordinates
        if center_y < 1 || center_y > img_height || center_x < 1 || center_x > img_width
            error('Kernel center (%d,%d) is outside image boundaries', center_y, center_x);
        end
        
        if interactive
            % Create figure for kernel size adjustment
            fig = figure('Name', sprintf('Kernel %d - Drag to resize, press Enter to confirm', n), ...
                        'Position', [100 100 500 500]);
            
            % Display the image
            imagesc(Y);
            colormap(gray);
            axis square;
            hold on;
            
            % Draw initial rectangle with target size
            half_height = floor(target_kernel_size(n,1)/2);
            half_width = floor(target_kernel_size(n,2)/2);
            
            y1 = max(1, center_y - half_height);
            y2 = min(img_height, center_y + half_height);
            x1 = max(1, center_x - half_width);
            x2 = min(img_width, center_x + half_width);
            
            % Create rectangle with center fixed
            h_rect = imrect(gca, [x1 y1 x2-x1 y2-y1]);
            
            % Create center dot indicator
            h_dot = plot(center_x, center_y, 'r.', 'MarkerSize', 20);
            
            % Add listener to update rectangle size while keeping center fixed
            addNewPositionCallback(h_rect, @(pos) updateRectangleSize(pos, center_x, center_y));
            
            % Set position constraints to maintain square aspect ratio
            setPositionConstraintFcn(h_rect, @(pos) constrainToSquare(pos, center_x, center_y, img_width, img_height));
            
            % Wait for user to finish (press Enter)
            position = wait(h_rect);
            
            % Get final position and size
            kernel_sizes(n,:) = [position(4), position(3)];  % [height, width]
            
            % Extract region with final size
            y1 = round(position(2));
            y2 = round(position(2) + position(4));
            x1 = round(position(1));
            x2 = round(position(1) + position(3));
            
            % Delete the center dot and close figure
            delete(h_dot);
            close(fig);
        else
            % Use target kernel size without interactive adjustment
            kernel_sizes(n,:) = target_kernel_size(n,:);
            
            % Calculate region boundaries
            half_height = floor(target_kernel_size(n,1)/2);
            half_width = floor(target_kernel_size(n,2)/2);
            
            y1 = max(1, center_y - half_height);
            y2 = min(img_height, center_y + half_height);
            x1 = max(1, center_x - half_width);
            x2 = min(img_width, center_x + half_width);
        end
        
        % Extract region
        selected_kernel = Y(y1:y2, x1:x2);
        
        % Project onto the oblique manifold and apply window
        A1{n} = proj2oblique(selected_kernel);
        disp(size(A1{n}));
        A1{n} = apply_window(A1{n}, window_type);
        
        fprintf('Kernel %d initialized at center (%d,%d) with size [%d,%d]\n', ...
            n, center_y, center_x, kernel_sizes(n,1), kernel_sizes(n,2));
    end
end

function new_pos = updateRectangleSize(pos, center_x, center_y)
    % Update rectangle size while keeping center fixed
    new_pos = pos;
    size = max(pos(3), pos(4));  % Use the larger of width or height
    new_pos(1) = center_x - size/2;
    new_pos(2) = center_y - size/2;
    new_pos(3) = size;
    new_pos(4) = size;
end

function new_pos = constrainToSquare(pos, center_x, center_y, img_width, img_height)
    % Constrain the rectangle to maintain square aspect ratio and stay within image bounds
    size = max(pos(3), pos(4));  % Use the larger of width or height
    
    % Calculate maximum possible size based on distance to image boundaries
    max_size_x = 2 * min(center_x - 1, img_width - center_x);
    max_size_y = 2 * min(center_y - 1, img_height - center_y);
    max_size = min(max_size_x, max_size_y);
    
    % Ensure size doesn't exceed image boundaries
    size = min(size, max_size);
    
    % Calculate new position keeping center fixed
    new_pos = [center_x - size/2, center_y - size/2, size, size];
end

function kernel_out = apply_window(kernel, window_type)
    % Helper function to apply window
    if iscell(window_type)
        % For windows with parameters (gaussian, kaiser)
        kernel_out = windowToKernel(kernel, window_type{1}, window_type{2});
    else
        % For windows without parameters
        kernel_out = windowToKernel(kernel, window_type);
    end
end