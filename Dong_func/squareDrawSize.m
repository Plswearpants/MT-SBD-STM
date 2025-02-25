function [square_size, position, mask] = squareDrawSize(matrix)
    % Display the matrix as an image
    figure;
    imshow(matrix, []);
    title('Determine the kernel size, drag a square covering the whole scattering pattern centered at a defect. Double-click to confirm');
    
    
    % Create the interactive rectangle ROI
    h = drawrectangle('DrawingArea', 'auto', 'FixedAspectRatio', true);
    
    % Add a listener to execute a callback function when the rectangle is finished drawing
    addlistener(h, 'ROIClicked', @(src, evt) clickCallback(src, evt, h));

    % Wait for user to double-click inside the ROI to confirm the selection
    uiwait;
    
    % Retrieve the final position of the ROI
    position = h.Position;
    
    % Calculate the side length of the square
    side_length = position(3);  % Since it's a square, width equals height
    
    % Draw the final square
    hold on;
    rectangle('Position', position, 'EdgeColor', 'r', 'LineWidth', 2);
    hold off;
    
    % Calculate and output the size of the square (side length) in pixels
    square_size = [round(side_length),round(side_length)];  % Ensure the size is in whole pixels
    disp(['Size of the square (side length in pixels): ', num2str(square_size)]);
    
    % Create the mask
    mask = create_mask(matrix, position);
    
    % Close the figure
    close;
end

function clickCallback(~, evt, h)
    % Resume execution only if the mouse click type is 'double'
    if strcmp(evt.SelectionType, 'double')
        uiresume;
    end
end

function mask = create_mask(matrix, position)
    % Create a mask with the same dimensions as the input matrix
    mask = zeros(size(matrix));
    
    % Get the coordinates and size of the rectangle
    x = round(position(1));
    y = round(position(2));
    width = round(position(3));
    height = round(position(4));
    
    % Set the region inside the rectangle to 1
    mask(y:y+height-1, x:x+width-1) = 1;
end
