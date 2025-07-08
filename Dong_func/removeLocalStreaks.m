function [corrected_data, streak_mask] = removeLocalStreaks(Y, slice_idx, min_value, max_streak_width)
%REMOVELOCALSTREAKS Interactive streak removal tool using multi-scale Laplacian detection
%   [corrected_data, streak_mask] = removeLocalStreaks(Y, slice_idx, min_value, max_streak_width) processes the specified slice of 3D data Y
%   using multi-scale Laplacian-based detection and neighbor interpolation correction
%
%   Inputs:
%       Y - 3D data array
%       slice_idx - Index of the slice to process (default: 150)
%       min_value - Optional minimum value for streak detection (default: interactive)
%       max_streak_width - Maximum streak width to detect (default: 3)
%
%   Outputs:
%       corrected_data - The corrected image data
%       streak_mask - Binary mask indicating detected streaks
%
%   Example:
%       [corrected, mask] = removeLocalStreaks(Y, 150);
%       [corrected, mask] = removeLocalStreaks(Y, 150, 0.5, 5);

% Initialize
if nargin < 2
    slice_idx = 150;
end
if nargin < 3
    min_value = [];
end
if nargin < 4
    max_streak_width = 3;
end
corrected_data = [];
streak_mask = [];

% Process data
data = Y(:,:,slice_idx);
[rows, cols] = size(data);

% Compute global mean for fallback
global_mean = mean(data(:));

% Compute multi-scale Laplacians
L_combined = zeros(size(data));
L_mag_combined = zeros(size(data));

for n = 1:max_streak_width
    % Shift left and right by n pixels
    data_left_n = [zeros(rows, n), data(:, 1:end-n)];
    data_right_n = [data(:, n+1:end), zeros(rows, n)];
    
    % Compute Laplacian at distance n
    L_n = data_left_n + data_right_n - 2*data;
    L_mag_n = abs(L_n);
    
    % Combine with previous scales (take maximum magnitude)
    L_mag_combined = max(L_mag_combined, L_mag_n);
end

% Create figure and store its handle
h_fig = figure('Name', 'Multi-Scale X-Direction Laplacian Streak Removal Analysis', 'Position', [100, 100, 1200, 800]);

% Plot original data
subplot(2, 2, 1);
imagesc(data);
title('Original Data');
axis square;
colormap gray;
colorbar;

% Plot binary mask (initially all zeros)
subplot(2, 2, 2);
h_mask = imagesc(zeros(size(L_mag_combined)));
title('Detected Streak Mask');
axis square;
colormap gray;
colorbar;

% Plot histogram with dual y-axis
subplot(2, 2, 3);
yyaxis left
h_hist = histogram(L_mag_combined);
title('Histogram of Combined Laplacian Magnitude');
axis square;
hold on;
h_min_line = xline(min(L_mag_combined(:)), 'r-', 'LineWidth', 2);
h_max_line = xline(max(L_mag_combined(:)), 'r-', 'LineWidth', 2);
ylabel('Count');

yyaxis right
h_var_line = plot([min(L_mag_combined(:)), max(L_mag_combined(:))], [0, 0], 'b-', 'LineWidth', 2);
ylabel('Variance');
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

% Create all UI controls first
min_label = uicontrol(panel, 'Style', 'text', 'String', 'Min:', 'Position', [10, 5, 40, 20]);

% Initialize slider value based on whether min_value was provided
initial_min = min(L_mag_combined(:));
if ~isempty(min_value)
    initial_min = min_value;
end

min_slider = uicontrol(panel, 'Style', 'slider', ...
    'Min', min(L_mag_combined(:)), ...
    'Max', max(L_mag_combined(:)), ...
    'Value', initial_min, ...
    'Position', [60, 5, 200, 20]);

min_edit = uicontrol(panel, 'Style', 'edit', ...
    'String', sprintf('%.3f', initial_min), ...
    'Position', [270, 5, 80, 20]);

done_button = uicontrol(panel, 'Style', 'pushbutton', ...
    'String', 'Done', ...
    'Position', [450, 5, 100, 40]);

% Set up callbacks after all controls are created
set(min_slider, 'Callback', @(src,event) updateContrast(src, event, h_mask, min_edit, h_corrected, data, L_mag_combined, h_min_line, h_max_line, done_button, h_var_line, global_mean));
set(min_edit, 'Callback', @(src,event) updateFromText(src, event, min_slider, h_mask, h_corrected, data, L_mag_combined, h_min_line, h_max_line, done_button, h_var_line, global_mean));
set(done_button, 'Callback', @(src,event) finish(src, h_corrected, h_mask));

% Initialize display
updateContrast(min_slider, [], h_mask, min_edit, h_corrected, data, L_mag_combined, h_min_line, h_max_line, done_button, h_var_line, global_mean);

% Wait for the figure to be closed
waitfor(h_fig);

% Get results from the base workspace
if evalin('base', 'exist(''temp_corrected_data'', ''var'')')
    corrected_data = evalin('base', 'temp_corrected_data');
    streak_mask = evalin('base', 'temp_streak_mask');
    evalin('base', 'clear temp_corrected_data temp_streak_mask');
end

end

function updateFromText(src, ~, min_slider, h_mask, h_corrected, data, L_mag_combined, h_min_line, h_max_line, done_button, h_var_line, global_mean)
    % Get the value from the text input
    new_val = str2double(get(src, 'String'));
    
    % Validate the input
    if isnan(new_val)
        % If invalid input, reset to current slider value
        set(src, 'String', sprintf('%.3f', get(min_slider, 'Value')));
        return;
    end
    
    % Clamp the value to the valid range
    new_val = max(min(L_mag_combined(:)), min(max(L_mag_combined(:)), new_val));
    
    % Update the slider
    set(min_slider, 'Value', new_val);
    
    % Update the display
    updateContrast(min_slider, [], h_mask, src, h_corrected, data, L_mag_combined, h_min_line, h_max_line, done_button, h_var_line, global_mean);
    
    % Force immediate update of all figures
    drawnow;
    
    % Update the text field with the clamped value if it was changed
    if new_val ~= str2double(get(src, 'String'))
        set(src, 'String', sprintf('%.3f', new_val));
    end
end

function finish(src, h_corrected, h_mask)
    % Get the current values from the button's UserData
    user_data = get(src, 'UserData');
    min_val = user_data(1).min_val;
    max_val = user_data(1).max_val;
    
    % Get the results
    corrected_data = get(h_corrected, 'CData');
    streak_mask = get(h_mask, 'CData') >= min_val & get(h_mask, 'CData') <= max_val;
    
    % Store results in base workspace
    assignin('base', 'temp_corrected_data', corrected_data);
    assignin('base', 'temp_streak_mask', streak_mask);
    
    % Close the figure
    close(gcf);
end

function updateContrast(src, ~, h_mask, min_edit, h_corrected, data, L_mag_combined, h_min_line, h_max_line, done_button, h_var_line, global_mean)
    % Get current values
    min_val = get(src, 'Value');
    max_val = max(L_mag_combined(:));
    
    % Update display
    caxis(h_mask.Parent, [min_val, max_val]);
    set(min_edit, 'String', sprintf('%.3f', min_val));
    set(h_min_line, 'Value', min_val);
    set(h_max_line, 'Value', max_val);
    
    % Find streaks
    streak_mask = L_mag_combined >= min_val & L_mag_combined <= max_val;
    [streak_rows, streak_cols] = find(streak_mask);
    
    % Correct streaks using left-right interpolation
    corrected = data;
    unique_cols = unique(streak_cols);
    
    for col = unique_cols'
        rows_in_col = streak_rows(streak_cols == col);
        
        % Process all streak points in this column
        if ~isempty(rows_in_col)
            left_col = col - 1;
            right_col = col + 1;
            left_vals = [];
            right_vals = [];
            if left_col >= 1
                left_vals = mean(data(rows_in_col, left_col));
            end
            if right_col <= size(data, 2)
                right_vals = mean(data(rows_in_col, right_col));
            end
            streak_avg = mean(data(rows_in_col, col));
            % Calculate expected value
            if ~isempty(left_vals) && ~isempty(right_vals)
                % If both neighbors exist, do linear interpolation
                left_avg = left_vals;
                right_avg = right_vals;
                expected_avg = left_avg + (col-left_col)*(right_avg-left_avg)/(right_col-left_col);
            elseif ~isempty(left_vals)
                expected_avg = left_vals;
            elseif ~isempty(right_vals)
                expected_avg = right_vals;
            else
                expected_avg = global_mean;
            end
            corrected(rows_in_col, col) = data(rows_in_col, col) - streak_avg + expected_avg;
        end
    end
    
    % Calculate variance of corrected data
    var_corrected = var(corrected(:));
    
    % Update variance line
    set(h_var_line, 'YData', [var_corrected, var_corrected]);
    
    % Update display
    set(h_corrected, 'CData', corrected);
    caxis(h_corrected.Parent, [min(data(:)), max(data(:))]);
    set(done_button, 'UserData', struct('min_val', min_val, 'max_val', max_val));
    
    % Force immediate update
    drawnow;
    
    % Update binary mask visualization
    set(h_mask, 'CData', streak_mask);
end


