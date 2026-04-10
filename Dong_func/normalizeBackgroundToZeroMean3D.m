function [normalized_data_3d, background_info, comment, position_used] = normalizeBackgroundToZeroMean3D(data_3d, rangeType, selected_slice, position)
% Normalizes each slice by subtracting the mean of a background area (selected once, reused).
%   If position is provided, that rectangle is used for all slices (no UI).
%   Otherwise the user selects a rectangle on selected_slice; it is returned for reuse.
%
% Arguments:
%   data_3d          3D array containing the data to be processed.
%   rangeType        Type of range for visualization ('global' or 'dynamic').
%   selected_slice   Slice index used only when selecting the region (default 1).
%   position         (Optional) [x, y, width, height] from a previous call. If provided, no UI.
%
% Returns:
%   normalized_data_3d   The 3D data with each slice normalized by the background mean.
%   background_info      [mean, variance] of the selected background area (last slice).
%   comment               Comment for logging.
%   position_used         The rectangle [x, y, width, height] used (for storing/replay).

arguments
    data_3d
    rangeType
    selected_slice = 1
    position = []  % [x, y, width, height]; if empty, user selects interactively
end

[rows, cols, k] = size(data_3d);

if isempty(position)
    figure;
    imagesc(data_3d(:,:,selected_slice));
    axis square;
    title('Select Background Area (This will be applied to all slices)');
    h = imrect;
    position = wait(h);
    close;
end

position_used = position;

x_min = round(position(1));
y_min = round(position(2));
x_max = round(position(1) + position(3));
y_max = round(position(2) + position(4));

% Clamp to image bounds
x_min = max(1, min(x_min, cols));
x_max = max(1, min(x_max, cols));
y_min = max(1, min(y_min, rows));
y_max = max(1, min(y_max, rows));
if x_min > x_max, [x_min, x_max] = deal(x_max, x_min); end
if y_min > y_max, [y_min, y_max] = deal(y_max, y_min); end

comment = sprintf("normalizeBackgroundToZeroMean3D(selected_slice:%d, region:[%d %d %d %d])|", ...
    selected_slice, x_min, y_min, x_max, y_max);

normalized_data_3d = zeros(size(data_3d));
background_info = [0, 0];

for i = 1:k
    background = data_3d(y_min:y_max, x_min:x_max, i);
    background_mean = mean(background(:));
    background_info = [background_mean, var(background(:))];
    data = data_3d(:, :, i);
    normalized_data_3d(:, :, i) = data - background_mean;
end
end
