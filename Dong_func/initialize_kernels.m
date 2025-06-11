function A1 = initialize_kernels(Y, num_kernels, kernel_sizes, kerneltype, window_type)
    % Inputs:
    %   Y: input image
    %   num_kernels: number of kernels to initialize
    %   kernel_size: size of each kernel [h w]
    %   kerneltype: 'random' or 'selected'
    %   window_type: (optional) type of window to apply
    %               empty or '' for no window
    %               'hann', 'hamming', 'blackman', 'gaussian', 'kaiser'
    %               For gaussian/kaiser, use cell array: {'gaussian', alpha} or {'kaiser', beta}
    % Outputs:
    %   A1: cell array of initialized kernels
    %       Each cell contains a kernel of size [h w] projected onto oblique manifold
    %       If window_type is specified, window function is applied to each kernel
    % Set default window_type to empty if not provided
    if nargin < 5
        window_type = '';
    end
    
    A1 = cell(1, num_kernels);

    switch lower(kerneltype)
        case 'random'
            for n = 1:num_kernels
                A1{n} = proj2oblique(randn(kernel_sizes(n,:)));
                % Apply window if specified
                if ~isempty(window_type)
                    A1{n} = apply_window(A1{n}, window_type);
                end
            end

        case 'selected'
            % Create figure for selection
            fig = figure;
            imagesc(Y);
            colorbar;
            axis square;
            title('Select initial kernel regions (Center dot indicates kernel center)');
            hold on;

            [img_height, img_width] = size(Y);

            for n = 1:num_kernels
                % Initial position
                init_x = 1;
                init_y = 1;
                rect_width = min(kernel_sizes(n,2), img_width);
                rect_height = min(kernel_sizes(n,1), img_height);
                
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
                    A1{n} = apply_window(A1{n}, window_type);
                end
                
                % Delete the center dot
                delete(h_dot);
            end
            close(fig);

        otherwise
            error('Invalid kernel initialization type. Choose "random" or "selected".');
    end
end

function new_pos = constrainPosition(pos, img_width, img_height, kernel_size)
    % Constrain the position of the rectangle to stay within the image bounds
    new_pos = pos;
    new_pos(1) = max(1, min(pos(1), img_width - kernel_size(2) + 1));
    new_pos(2) = max(1, min(pos(2), img_height - kernel_size(1) + 1));
    new_pos(3) = kernel_size(2);
    new_pos(4) = kernel_size(1);
end

function updateCenterDot(pos, h_dot)
    % Update center dot position when rectangle moves
    center_x = pos(1) + pos(3)/2;
    center_y = pos(2) + pos(4)/2;
    set(h_dot, 'XData', center_x, 'YData', center_y);
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
