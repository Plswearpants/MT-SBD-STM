function [mask] = linesegmentEdge2Edge(data, direction, varargin)
% LINESEGMENTEDGE2EDGE Draw line segments from edge to edge with user input and create a mask
%   mask = linesegmentEdge2Edge(data, direction)
%   mask = linesegmentEdge2Edge(data, direction, 'LineWidth', width)
%
% Inputs:
%   data: 2D or 3D data array (if 3D, first slice is used for display)
%   direction: 'horizontal' or 'vertical'
%   'LineWidth': (optional) Width of the line segments in pixels (default: 1)
%
% Outputs:
%   mask: Binary mask where 1 indicates points within the line segments

% Parse optional inputs
p = inputParser;
addParameter(p, 'LineWidth', 1, @(x) isnumeric(x) && x > 0);
parse(p, varargin{:});
line_width = p.Results.LineWidth;

% if data is 3D, use d3gridDisplay to display and let user choose the slice
if ndims(data) == 3
    d3gridDisplay(data, 'dynamic');
end

% prompt user to choose the slice
slice = input('Please choose the slice number: ');
close;

% Display the slice
figure;
if ndims(data) == 3
    imagesc(data(:,:,slice));
else
    imagesc(data);
end
title('Click to add points. Press Enter when done.');
colormap(gray);
axis equal;
hold on;

% Get image dimensions
[rows, cols] = size(data(:,:,1));

% Initialize arrays to store points
x_points = [];
y_points = [];

% Get first point from user
[x0, y0] = ginput(1);

% Add 0th point (edge point) based on direction
switch direction
    case 'horizontal'
        x_points = [0; x0];  % Start from left edge
        y_points = [y0; y0]; % Same y as first clicked point
    case 'vertical'
        x_points = [x0; x0]; % Same x as first clicked point
        y_points = [0; y0];  % Start from top edge
end

% Plot points and line
plot(x_points(1), y_points(1), 'ro', 'MarkerSize', 8);
plot(x_points(2), y_points(2), 'ro', 'MarkerSize', 8);
plot(x_points, y_points, 'r-', 'LineWidth', 2);

% Let user add intermediate points
while true
    % Get next point
    [x, y] = ginput(1);
    
    % Check if user pressed Enter (empty input)
    if isempty(x)
        break;
    end
    
    % Add point to arrays
    x_points = [x_points; x];
    y_points = [y_points; y];
    
    % Plot point and line to previous point
    plot(x, y, 'ro', 'MarkerSize', 8);
    plot([x_points(end-1), x_points(end)], [y_points(end-1), y_points(end)], 'r-', 'LineWidth', 2);
end

% Add last point based on direction
switch direction
    case 'horizontal'
        x_last = cols;
        y_last = y_points(end);  % Use same y as last clicked point
    case 'vertical'
        x_last = x_points(end);  % Use same x as last clicked point
        y_last = rows;
end

% Add last point to arrays
x_points = [x_points; x_last];
y_points = [y_points; y_last];

% Plot last point and final line segment
plot(x_last, y_last, 'ro', 'MarkerSize', 8);
plot([x_points(end-1), x_points(end)], [y_points(end-1), y_points(end)], 'r-', 'LineWidth', 2);

% Create binary mask
mask = false(rows, cols);

% Create coordinate grids
[X, Y] = meshgrid(1:cols, 1:rows);

% For each line segment
for i = 1:length(x_points)-1
    % Get points along the line segment
    [xq, yq] = get_points_along_line([x_points(i), x_points(i+1)], [y_points(i), y_points(i+1)]);
    
    % For each point along the line segment
    for j = 1:length(xq)
        % Calculate distances to all grid points
        distances = sqrt((X - xq(j)).^2 + (Y - yq(j)).^2);
        
        % Set mask to true for points within the specified width of the line
        mask(distances <= line_width/2) = true;
    end
end

% Display the mask
figure;
imagesc(mask);
title(sprintf('Generated Mask (Line Width: %d pixels)', line_width));
colormap(gray);
axis equal;

end

function [xq, yq] = get_points_along_line(x, y)
    % Get points along a line segment with sub-pixel resolution
    num_points = max(abs(diff(x)), abs(diff(y))) * 2;
    t = linspace(0, 1, num_points);
    xq = interp1([0 1], x, t);
    yq = interp1([0 1], y, t);
end 