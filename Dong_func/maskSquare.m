function [mask] = maskSquare(Y, center, slice, shape)
% maskSquare allows user to draw a square or rectangular ROI on 2D/3D data
% Input:
%   Y - 2D or 3D data matrix
%   center - (optional) logical, if true enforces same distance to edges
%            between the square and the original image (default: false)
%   slice - (optional) for 3D data, which slice to display
%   shape - (optional) 'square' (default) or 'rectangle'
% Output:
%   mask - binary mask of the selected region

% Set default values for center, slice, and shape if not provided
if nargin < 2
    center = false;
end
if nargin < 3
    slice = 1;
end
if nargin < 4
    shape = 'square';
end

% Display the first slice of the data
figure;
if ndims(Y) == 3
    imagesc(Y(:,:,slice));
else
    imagesc(Y);
end
axis image;
%colormap hot;
title('Draw a ROI');

% Let user draw a square or rectangular ROI
if strcmpi(shape, 'rectangle')
    h = drawrectangle('Color', 'g', 'FixedAspectRatio', false);
else
    h = drawrectangle('Color', 'r', 'FixedAspectRatio', true);
end
wait(h);

% Get ROI properties
pos = h.Position;  % [x y width height]
x = round(pos(1));
y = round(pos(2));
width = round(pos(3));
height = round(pos(4));

% If center is true and shape is square, adjust position to maintain equal distance to edges
if center && strcmpi(shape, 'square')
    [rows, cols] = size(Y(:,:,1));
    % Calculate current distances to edges
    dist_left = x;
    dist_right = cols - (x + width - 1);
    dist_top = y;
    dist_bottom = rows - (y + height - 1);
    % Calculate average distance to maintain
    avg_dist_x = round((dist_left + dist_right) / 2);
    avg_dist_y = round((dist_top + dist_bottom) / 2);
    % Adjust position to maintain equal distance
    x = avg_dist_x;
    y = avg_dist_y;
    % Ensure the square fits within the image bounds
    if x + width > cols
        width = cols - x;
    end
    if y + height > rows
        height = rows - y;
    end
end

% Create mask with proper rectangular format
mask = false(size(Y,1), size(Y,2));
mask(y:y+height-1, x:x+width-1) = true;

end 