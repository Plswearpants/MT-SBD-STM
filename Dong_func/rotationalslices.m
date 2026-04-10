function [data2D, angles, radius] = rotationalslices(data3D, rangetype, line_width, radius, naked, energy_range)
% ROTATIONAL_SLICES Creates 2D slices from 3D data by rotating a line through the center
%   This function creates a series of 2D slices by rotating a line through the center
%   of the data at different angles. The user adjusts a pre-drawn circle to determine the line length.
%
% Inputs:
%   data3D - 3D dataset to visualize
%   rangetype - 'dynamic' (normalized per slice), 'global' (raw data), or 'log' (logarithmic scaling)
%   line_width - Width of the line in pixels (default: 1)
%   radius - Optional radius for the circle. If not provided, user will be prompted to draw a circle
%   naked - If true, hide plot decorations but keep interactive controls
%   energy_range - Optional [Emin, Emax] in meV for slice-energy readout
%
% Outputs:
%   data2D - 3D array where:
%            dim1: position along the line
%            dim2: original data slices (energy)
%            dim3: angle slices
%   angles - Array of angles used for rotation (in radians)
%   radius - Radius of the circle used for slicing
%
% Example:
%   [data2D, angles] = rotational_slices(LDoS_noisy, 'log', 3)
%   [data2D, angles, radius] = rotational_slices(LDoS_noisy, 'log', 3, 50)
%
% Created: March 2024

% Validate inputs
arguments
    data3D {mustBeNumeric}
    rangetype {mustBeMember(rangetype, {'dynamic', 'global', 'log'})} = 'dynamic'
    line_width {mustBeInteger, mustBePositive} = 1
    radius {mustBeNumeric} = []
    naked (1,1) logical = false
    energy_range (1,2) double = [NaN NaN]
end

% Calculate center point
[rows, cols, ~] = size(data3D);
center = [floor(rows/2)+1, floor(cols/2)+1];
num_slices = size(data3D,3);

% Build energy axis for readout; fallback to slice index when not provided.
if all(isfinite(energy_range))
    energy_axis = linspace(energy_range(1), energy_range(2), num_slices);
    energy_unit = 'meV';
else
    energy_axis = 1:num_slices;
    energy_unit = 'idx';
end

% If radius not provided, get it interactively
if isempty(radius)
    % Display the first slice
    figure;
    imagesc(data3D(:,:,floor(size(data3D,3)/2)));
    colormap gray;
    axis equal;
    title('Adjust circle radius to determine line length');
    hold on;

    % Create initial circle with radius 1
    h = drawcircle('Center', [center(2), center(1)], 'Radius', 1);
    h.FaceAlpha = 0.1;  % Make circle semi-transparent

    % Wait for user to finish adjusting
    wait(h);

    % Get circle properties
    radius = h.Radius;
    
    % Clean up
    delete(h);
    hold off;
end

% Generate angles (0 to pi)
angles = linspace(0, pi, 360);

% Calculate maximum possible points for any angle
% This will be the diameter of the circle in pixels
max_points = 2 * round(radius);

% Pre-allocate output array with correct dimensions
data2D = zeros(max_points, size(data3D, 3), length(angles));

% Process each angle
for i = 1:length(angles)
    % Calculate endpoints on circle
    theta = angles(i);
    pointA = round([center(1) + radius*cos(theta), center(2) + radius*sin(theta)]);
    pointB = round([center(1) - radius*cos(theta), center(2) - radius*sin(theta)]);
    
    % Ensure points are within bounds
    pointA = max(1, min([rows, cols], pointA));
    pointB = max(1, min([rows, cols], pointB));
    
    % Calculate perpendicular direction for line width
    line_vec = [pointB(2) - pointA(2), pointB(1) - pointA(1)];
    line_vec = line_vec / norm(line_vec);
    perp_vec = [-line_vec(2), line_vec(1)];  % Rotate 90 degrees
    
    % Initialize array for width-averaged data
    width_data = zeros(max_points, size(data3D, 3), line_width);
    
    % Process each parallel line
    for w = 1:line_width
        % Calculate offset for this line
        offset = (w - (line_width+1)/2) * perp_vec;
        
        % Calculate offset points
        pointA_offset = round(pointA + offset);
        pointB_offset = round(pointB + offset);
        
        % Ensure points are within bounds
        pointA_offset = max(1, min([rows, cols], pointA_offset));
        pointB_offset = max(1, min([rows, cols], pointB_offset));
        
        % Get mask using gridMaskLineSegment
        [mask, ~, ~, ~, ~] = gridMaskLineSegment(data3D(:,:,1), pointA_offset, pointB_offset);
        
        if any(mask(:))
            % Get coordinates of masked points
            [y_coords, x_coords] = find(mask);
            
            % Calculate distances along the line
            distances = sqrt((x_coords - pointA_offset(2)).^2 + (y_coords - pointA_offset(1)).^2);
            [~, sort_idx] = sort(distances);
            x_coords = x_coords(sort_idx);
            y_coords = y_coords(sort_idx);
            
            % Create 2D data array for current angle
            switch rangetype
                case 'dynamic'
                    for z = 1:size(data3D, 3)
                        values = data3D(sub2ind(size(data3D), y_coords, x_coords, z*ones(size(x_coords))));
                        % Interpolate to fixed number of points
                        interpolated = interp1(1:length(values), values, ...
                                          linspace(1,length(values),max_points));
                        % Normalize to range 0-1
                        width_data(:,z,w) = (interpolated - min(interpolated)) / (max(interpolated) - min(interpolated));
                    end
                case 'global'
                    for z = 1:size(data3D, 3)
                        values = data3D(sub2ind(size(data3D), y_coords, x_coords, z*ones(size(x_coords))));
                        % Interpolate to fixed number of points
                        interpolated = interp1(1:length(values), values, ...
                                          linspace(1,length(values),max_points));
                        % no normalization
                        width_data(:,z,w) = interpolated;
                    end
                case 'log'
                    for z = 1:size(data3D, 3)
                        values = data3D(sub2ind(size(data3D), y_coords, x_coords, z*ones(size(x_coords))));
                        % Interpolate to fixed number of points
                        interpolated = interp1(1:length(values), values, ...
                                          linspace(1,length(values),max_points));
                        % Apply logarithmic scaling
                        % Add small offset to avoid log(0)
                        offset = 1e-10;
                        width_data(:,z,w) = log(abs(interpolated));
                    end
            end
        end
    end
    
    % Average the data from all parallel lines
    data2D(:,:,i) = mean(width_data, 3);
end

% Create figure with two subplots
fig = figure('Position', [100 100 1200 560]);

% Initial slice index for embedded 3D viewer behavior
slice_idx = 1;

% First subplot: Original data with rotating line
subplot(1,2,1);
h1 = imagesc(data3D(:,:,slice_idx));
colormap invgray;
axis equal;
hold on;
h_line = plot([pointA(2), pointB(2)], [pointA(1), pointB(1)], 'r-', 'LineWidth', 2);
if ~naked
    title(sprintf('Original data with line (slice %d/%d)', slice_idx, size(data3D,3)));
else
    axis off;
end

% Second subplot: Rotational slice data
subplot(1,2,2);
h2 = imagesc(data2D(:,:,1)');
colormap invgray;
set(gca, 'YDir', 'normal'); % Make slice 1 appear at bottom
hold on;
h_slice_line = yline(slice_idx, 'r-', 'LineWidth', 1.5);
hold off;
if ~naked
    title(sprintf('Angle: %.2f° (0° = vertical)', rad2deg(angles(1))));
    xlabel('Position along line');
    ylabel('Energy slice');
    colorbar;
else
    axis off;
end

% Store data in figure for callback
setappdata(fig, 'data2D', data2D);
setappdata(fig, 'data3D', data3D);
setappdata(fig, 'angles', angles);
setappdata(fig, 'center', center);
setappdata(fig, 'radius', radius);
setappdata(fig, 'h1', h1);
setappdata(fig, 'h2', h2);
setappdata(fig, 'h_line', h_line);
setappdata(fig, 'h_slice_line', h_slice_line);
setappdata(fig, 'line_width', line_width);
setappdata(fig, 'slice_idx', slice_idx);
setappdata(fig, 'naked', naked);
setappdata(fig, 'energy_axis', energy_axis);
setappdata(fig, 'energy_unit', energy_unit);

% Keep basic slider functionality available in both normal and naked modes
slider = uicontrol('Style', 'slider',...
    'Min', 1, 'Max', length(angles),...
    'Value', 1,...
    'Position', [20 20 400 20],...
    'Callback', @updatePlot);

slice_slider = uicontrol('Style', 'slider',...
    'Min', 1, 'Max', size(data3D,3),...
    'Value', slice_idx,...
    'SliderStep', [1/max(1,size(data3D,3)-1), min(10/max(1,size(data3D,3)-1),1)],...
    'Position', [560 20 400 20],...
    'Callback', @updateSlice);
setappdata(fig, 'slice_slider', slice_slider);

% Keep angle indicator visible in both normal and naked modes
angle_text = uicontrol('Style', 'text',...
    'Position', [430 20 100 20],...
    'String', sprintf('%.2f°', rad2deg(angles(1))));
% Keep slice indicator visible in both normal and naked modes
slice_text = uicontrol('Style', 'text',...
    'Position', [970 20 160 20],...
    'String', sprintf('Slice %d/%d', slice_idx, size(data3D,3)));
energy_text = uicontrol('Style', 'text',...
    'Position', [970 0 160 20],...
    'String', sprintf('Energy: %.1f %s', energy_axis(slice_idx), energy_unit));
line_toggle = uicontrol('Style', 'togglebutton',...
    'Position', [970 45 160 24],...
    'String', 'Lines: On',...
    'Value', 1,...
    'Callback', @toggleLines);
setappdata(fig, 'angle_text', angle_text);
setappdata(fig, 'slice_text', slice_text);
setappdata(fig, 'energy_text', energy_text);
setappdata(fig, 'line_toggle', line_toggle);

end

% Callback function to update plot
function updatePlot(src, ~)
    % Get data from figure
    fig = ancestor(src, 'figure');
    data2D = getappdata(fig, 'data2D');
    data3D = getappdata(fig, 'data3D');
    angles = getappdata(fig, 'angles');
    center = getappdata(fig, 'center');
    radius = getappdata(fig, 'radius');
    h1 = getappdata(fig, 'h1');
    h2 = getappdata(fig, 'h2');
    h_line = getappdata(fig, 'h_line');
    angle_text = getappdata(fig, 'angle_text');
    line_width = getappdata(fig, 'line_width');
    naked = getappdata(fig, 'naked');
    
    % Get current angle index
    idx = round(src.Value);
    theta = angles(idx);
    
    % Update line endpoints
    pointA = round([center(1) + radius*cos(theta), center(2) + radius*sin(theta)]);
    pointB = round([center(1) - radius*cos(theta), center(2) - radius*sin(theta)]);
    
    % Calculate perpendicular direction for line width
    line_vec = [pointB(2) - pointA(2), pointB(1) - pointA(1)];
    line_vec = line_vec / norm(line_vec);
    perp_vec = [-line_vec(2), line_vec(1)];  % Rotate 90 degrees
    
    % Update line in first subplot (note: imagesc uses [x,y] for plot)
    set(h_line, 'XData', [pointA(2), pointB(2)], 'YData', [pointA(1), pointB(1)]);
    
    % Update second subplot
    h2.CData = data2D(:,:,idx)';
    if ~naked
        title(sprintf('Angle: %.2f° (0° = vertical)', rad2deg(angles(idx))));
    end
    
    % Update angle text
    if ~isempty(angle_text) && isgraphics(angle_text)
        set(angle_text, 'String', sprintf('%.2f°', rad2deg(angles(idx))));
    end
end 

% Callback function to update displayed slice and horizontal marker
function updateSlice(src, ~)
    fig = ancestor(src, 'figure');
    data3D = getappdata(fig, 'data3D');
    h1 = getappdata(fig, 'h1');
    h_slice_line = getappdata(fig, 'h_slice_line');
    h_line = getappdata(fig, 'h_line');
    slice_text = getappdata(fig, 'slice_text');
    energy_text = getappdata(fig, 'energy_text');
    energy_axis = getappdata(fig, 'energy_axis');
    energy_unit = getappdata(fig, 'energy_unit');
    naked = getappdata(fig, 'naked');

    idx = round(src.Value);
    idx = max(1, min(size(data3D,3), idx));
    setappdata(fig, 'slice_idx', idx);

    % Update subplot 1 image content (embedded 3D-slice browsing behavior)
    h1.CData = data3D(:,:,idx);

    % Keep overlay line and update title
    ax1 = ancestor(h1, 'axes');
    axes(ax1); %#ok<LAXES>
    if ~naked
        title(ax1, sprintf('Original data with line (slice %d/%d)', idx, size(data3D,3)));
    end
    uistack(h_line, 'top');

    % Update subplot 2 horizontal line to corresponding y-position
    if ~isempty(h_slice_line) && isgraphics(h_slice_line)
        h_slice_line.Value = idx;
    end

    % Update slice text readout
    if ~isempty(slice_text) && isgraphics(slice_text)
        set(slice_text, 'String', sprintf('Slice %d/%d', idx, size(data3D,3)));
    end
    if ~isempty(energy_text) && isgraphics(energy_text)
        set(energy_text, 'String', sprintf('Energy: %.1f %s', energy_axis(idx), energy_unit));
    end
end

% Callback to toggle visibility of red line indicators
function toggleLines(src, ~)
    fig = ancestor(src, 'figure');
    h_line = getappdata(fig, 'h_line');
    h_slice_line = getappdata(fig, 'h_slice_line');

    if src.Value == 1
        if isgraphics(h_line)
            h_line.Visible = 'on';
        end
        if isgraphics(h_slice_line)
            h_slice_line.Visible = 'on';
        end
        src.String = 'Lines: On';
    else
        if isgraphics(h_line)
            h_line.Visible = 'off';
        end
        if isgraphics(h_slice_line)
            h_slice_line.Visible = 'off';
        end
        src.String = 'Lines: Off';
    end
end