function [corrected_data, streak_mask, streak_indices] = interpolateLocalStreaks(Y, slice_idx, min_value, provided_streak_indices)
%INTERPOLATELOCALSTREAKS Interactive streak interpolating tool using Laplacian and neighbor interpolation
%   [corrected_data, streak_mask, streak_indices] = interpolateLocalStreaks(Y, slice_idx, min_value, provided_streak_indices)
%   processes the specified slice of 3D data Y using two methods: Laplacian-based detection followed by neighbor interpolation
%
%   Inputs:
%       Y - 3D data array
%       slice_idx - Index of the slice to process (default: 150)
%       min_value - Optional minimum value for streak detection (default: interactive)
%       provided_streak_indices - Optional Nx2 array of [row,col] indices for streak points
%
%   Outputs:
%       corrected_data - The corrected image data
%       streak_mask - Binary mask indicating detected streaks
%       streak_indices - Nx2 array of [row,col] indices for streak points
%
%   Example:
%       [corrected, mask, indices] = interpolateLocalStreaks(Y, 150);
%       [corrected, mask, indices] = interpolateLocalStreaks(Y, 150, 0.5);
%       [corrected, mask] = interpolateLocalStreaks(Y, 150, [], provided_indices);

% Initialize
if nargin < 2
    slice_idx = 150;
end
if nargin < 3
    min_value = [];
end
corrected_data = [];
streak_mask = [];
streak_indices = [];

% Process data
data = single(Y(:,:,slice_idx));  % Convert to single precision
[rows, cols] = size(data);

% If streak indices are provided, directly apply interpolation
if nargin >= 4 && ~isempty(provided_streak_indices)
    corrected_data = data;
    streak_mask = false(size(data));
    streak_indices = provided_streak_indices;
    
    % Vectorized processing of provided streak points
    valid_streaks = false(size(data));
    valid_streaks(sub2ind(size(data), provided_streak_indices(:,1), provided_streak_indices(:,2))) = true;
    
    % Check for valid streaks (both neighbors are streak points)
    valid_streaks(:,1) = false;  % Remove first column
    valid_streaks(:,end) = false;  % Remove last column
    valid_streaks = valid_streaks & circshift(valid_streaks, [0 -1]) & circshift(valid_streaks, [0 1]);
    
    % Apply interpolation
    corrected_data(valid_streaks) = (data(circshift(valid_streaks, [0 -1])) + data(circshift(valid_streaks, [0 1]))) / 2;
    streak_mask = valid_streaks;
    return;
end

% Compute Laplacian efficiently
L = zeros(size(data), 'single');
L(:,2:end-1) = data(:,1:end-2) + data(:,3:end) - 2*data(:,2:end-1);
L_mag = abs(L);

% Create figure and store its handle
h_fig = figure('Name', 'X-Direction Laplacian Streak Removal Analysis', 'Position', [100, 100, 1200, 800]);

% Plot original data
subplot(2, 2, 1);
h_orig = imagesc(data);
title('Original Data');
axis square;
colormap gray;
colorbar;

% Plot binary mask (initially all zeros)
subplot(2, 2, 2);
h_mask = imagesc(false(size(L_mag)));
title('Detected Streak Mask');
axis square;
colormap gray;
colorbar;

% Plot histogram
subplot(2, 2, 3);
h_hist = histogram(L_mag);
title('Histogram of Laplacian Magnitude');
axis square;
hold on;
h_min_line = xline(min(L_mag(:)), 'r-', 'LineWidth', 2);
h_max_line = xline(max(L_mag(:)), 'r-', 'LineWidth', 2);
xlim([min(L_mag(:)), max(L_mag(:))/2]);
hold off;

% Plot corrected image
subplot(2, 2, 4);
h_corrected = imagesc(data);
title('Corrected Image');
axis square;
colormap gray;
colorbar;

% Create controls
panel = uipanel('Position', [0.1, 0.02, 0.8, 0.05]);

% Initialize slider value
initial_min = min(L_mag(:));
if ~isempty(min_value)
    initial_min = min_value;
end

% Min slider
uicontrol(panel, 'Style', 'text', 'String', 'Min:', 'Position', [10, 5, 40, 20]);
min_slider = uicontrol(panel, 'Style', 'slider', ...
    'Min', min(L_mag(:)), ...
    'Max', max(L_mag(:))/2, ...
    'Value', initial_min, ...
    'Position', [60, 5, 300, 20]);

% Add value display
min_text = uicontrol(panel, 'Style', 'text', ...
    'String', sprintf('%.3f', initial_min), ...
    'Position', [370, 5, 60, 20]);

% Done button
done_button = uicontrol(panel, 'Style', 'pushbutton', ...
    'String', 'Done', ...
    'Position', [450, 5, 100, 40], ...
    'Callback', @(src,event) finish(src, h_corrected, h_mask));

% Set up callbacks with debounce
set(min_slider, 'Callback', @(src,event) debouncedUpdate(src, event, h_mask, min_text, h_corrected, data, L_mag, h_min_line, h_max_line, done_button));

% Initialize display
updateContrast(min_slider, [], h_mask, min_text, h_corrected, data, L_mag, h_min_line, h_max_line, done_button);

% Wait for the figure to be closed
waitfor(h_fig);

% Get results from the base workspace
if evalin('base', 'exist(''temp_corrected_data'', ''var'')')
    corrected_data = evalin('base', 'temp_corrected_data');
    streak_mask = evalin('base', 'temp_streak_mask');
    streak_indices = evalin('base', 'temp_streak_indices');
    evalin('base', 'clear temp_corrected_data temp_streak_mask temp_streak_indices');
end

end

function finish(src, h_corrected, h_mask)
    % Get the current values from the button's UserData
    user_data = get(src, 'UserData');
    min_val = user_data(1).min_val;
    max_val = user_data(1).max_val;
    
    % Get the results
    corrected_data = get(h_corrected, 'CData');
    streak_mask = get(h_mask, 'CData');
    [streak_rows, streak_cols] = find(streak_mask);
    streak_indices = [streak_rows, streak_cols];
    
    % Store results in base workspace
    assignin('base', 'temp_corrected_data', corrected_data);
    assignin('base', 'temp_streak_mask', streak_mask);
    assignin('base', 'temp_streak_indices', streak_indices);
    
    % Close the figure
    close(gcf);
end

function debouncedUpdate(src, event, h_mask, min_text, h_corrected, data, L_mag, h_min_line, h_max_line, done_button)
    persistent lastUpdate
    if isempty(lastUpdate)
        lastUpdate = tic;
    end
    
    % Only update if 0.1 seconds have passed since last update
    if toc(lastUpdate) > 0.1
        updateContrast(src, event, h_mask, min_text, h_corrected, data, L_mag, h_min_line, h_max_line, done_button);
        lastUpdate = tic;
    end
end

function updateContrast(src, ~, h_mask, min_text, h_corrected, data, L_mag, h_min_line, h_max_line, done_button)
    % Get current values
    min_val = get(src, 'Value');
    max_val = max(L_mag(:));
    
    % Update binary mask visualization
    streak_mask = L_mag >= min_val & L_mag <= max_val;
    set(h_mask, 'CData', streak_mask);
    caxis(h_corrected.Parent, [min(data(:)), max(data(:))]);
    set(min_text, 'String', sprintf('%.3f', min_val));
    set(h_min_line, 'Value', min_val);
    set(h_max_line, 'Value', max_val);
    
    % Process valid streaks
    valid_streaks = streak_mask;
    valid_streaks(:,1) = false;  % Remove first column
    valid_streaks(:,end) = false;  % Remove last column
    valid_streaks = valid_streaks & circshift(streak_mask, [0 -1]) & circshift(streak_mask, [0 1]);
    
    % Apply interpolation
    corrected = data;
    corrected(valid_streaks) = (data(circshift(valid_streaks, [0 -1])) + data(circshift(valid_streaks, [0 1]))) / 2;
    
    % Update display
    set(h_corrected, 'CData', corrected);
    
    % Store values for finish function
    set(done_button, 'UserData', struct('min_val', min_val, 'max_val', max_val));
end


