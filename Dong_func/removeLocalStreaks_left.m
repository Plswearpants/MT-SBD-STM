function [corrected_data, streak_mask] = removeLocalStreaks_left(Y, slice_idx, min_value)
%REMOVELOCALSTREAKS_LEFT Interactive streak removal tool using Laplacian and left neighbor interpolation
%   [corrected_data, streak_mask] = removeLocalStreaks_left(Y, slice_idx, min_value) processes the specified slice of 3D data Y
%   using two methods: Laplacian-based correction followed by left neighbor interpolation
%
%   Inputs:
%       Y - 3D data array
%       slice_idx - Index of the slice to process (default: 150)
%       min_value - Optional minimum value for streak detection (default: interactive)
%
%   Outputs:
%       corrected_data - The corrected image data
%       streak_mask - Binary mask indicating detected streaks
%
%   Example:
%       [corrected, mask] = removeLocalStreaks_left(Y, 150);
%       [corrected, mask] = removeLocalStreaks_left(Y, 150, 0.5);

% Initialize
if nargin < 2
    slice_idx = 150;
end
if nargin < 3
    min_value = [];
end
corrected_data = [];
streak_mask = [];

% Process data
data = Y(:,:,slice_idx);
[rows, cols] = size(data);

% Compute Laplacian
L = zeros(size(data));
% Shift left and right
data_left = [zeros(rows,1), data(:,1:end-1)];
data_right = [data(:,2:end), zeros(rows,1)];
% Compute Laplacian using matrix operations
L = data_left + data_right - 2*data;
L_mag = abs(L);

% Create figure and store its handle
h_fig = figure('Name', 'X-Direction Laplacian Streak Removal Analysis', 'Position', [100, 100, 1200, 800]);

% Plot original data
subplot(2, 2, 1);
imagesc(data);
title('Original Data');
axis square;
colormap parula;
colorbar;

% Plot Laplacian
subplot(2, 2, 2);
h_plot = imagesc(L_mag);
title('X-Direction Laplacian Magnitude');
axis square;
colormap parula;
colorbar;

% Plot histogram with dual y-axis
subplot(2, 2, 3);
yyaxis left
h_hist = histogram(L_mag);
title('Histogram of Laplacian Magnitude');
axis square;
hold on;
h_min_line = xline(min(L_mag(:)), 'r-', 'LineWidth', 2);
h_max_line = xline(max(L_mag(:)), 'r-', 'LineWidth', 2);
ylabel('Count');

yyaxis right
h_var_line = plot([min(L_mag(:)), max(L_mag(:))], [0, 0], 'b-', 'LineWidth', 2);
ylabel('Variance');
hold off;

% Plot corrected image
subplot(2, 2, 4);
h_corrected = imagesc(data);
title('Corrected Image');
axis square;
colormap parula;
colorbar;

% Create controls
panel = uipanel('Position', [0.1, 0.02, 0.8, 0.05]);

% Create all UI controls first
min_label = uicontrol(panel, 'Style', 'text', 'String', 'Min:', 'Position', [10, 5, 40, 20]);

% Initialize slider value based on whether min_value was provided
initial_min = min(L_mag(:));
if ~isempty(min_value)
    initial_min = min_value;
end

min_slider = uicontrol(panel, 'Style', 'slider', ...
    'Min', min(L_mag(:)), ...
    'Max', max(L_mag(:)), ...
    'Value', initial_min, ...
    'Position', [60, 5, 200, 20]);

min_edit = uicontrol(panel, 'Style', 'edit', ...
    'String', sprintf('%.3f', initial_min), ...
    'Position', [270, 5, 80, 20]);

done_button = uicontrol(panel, 'Style', 'pushbutton', ...
    'String', 'Done', ...
    'Position', [450, 5, 100, 40]);

% Set up callbacks after all controls are created
set(min_slider, 'Callback', @(src,event) updateContrast(src, event, h_plot, min_edit, h_corrected, data, L_mag, h_min_line, h_max_line, done_button, h_var_line));
set(min_edit, 'Callback', @(src,event) updateFromText(src, event, min_slider, h_plot, h_corrected, data, L_mag, h_min_line, h_max_line, done_button, h_var_line));
set(done_button, 'Callback', @(src,event) finish(src, h_corrected, h_plot));

% Initialize display
updateContrast(min_slider, [], h_plot, min_edit, h_corrected, data, L_mag, h_min_line, h_max_line, done_button, h_var_line);

% Wait for the figure to be closed
waitfor(h_fig);

% Get results from the base workspace
if evalin('base', 'exist(''temp_corrected_data'', ''var'')')
    corrected_data = evalin('base', 'temp_corrected_data');
    streak_mask = evalin('base', 'temp_streak_mask');
    evalin('base', 'clear temp_corrected_data temp_streak_mask');
end

end

function updateFromText(src, ~, min_slider, h_plot, h_corrected, data, L_mag, h_min_line, h_max_line, done_button, h_var_line)
    % Get the value from the text input
    new_val = str2double(get(src, 'String'));
    
    % Validate the input
    if isnan(new_val)
        % If invalid input, reset to current slider value
        set(src, 'String', sprintf('%.3f', get(min_slider, 'Value')));
        return;
    end
    
    % Clamp the value to the valid range
    new_val = max(min(L_mag(:)), min(max(L_mag(:)), new_val));
    
    % Update the slider
    set(min_slider, 'Value', new_val);
    
    % Update the display
    updateContrast(min_slider, [], h_plot, src, h_corrected, data, L_mag, h_min_line, h_max_line, done_button, h_var_line);
    
    % Force immediate update of all figures
    drawnow;
    
    % Update the text field with the clamped value if it was changed
    if new_val ~= str2double(get(src, 'String'))
        set(src, 'String', sprintf('%.3f', new_val));
    end
end

function finish(src, h_corrected, h_plot)
    % Get the current values from the button's UserData
    user_data = get(src, 'UserData');
    min_val = user_data(1).min_val;
    max_val = user_data(1).max_val;
    
    % Get the results
    corrected_data = get(h_corrected, 'CData');
    streak_mask = get(h_plot, 'CData') >= min_val & get(h_plot, 'CData') <= max_val;
    
    % Store results in base workspace
    assignin('base', 'temp_corrected_data', corrected_data);
    assignin('base', 'temp_streak_mask', streak_mask);
    
    % Close the figure
    close(gcf);
end

function updateContrast(src, ~, h_plot, min_edit, h_corrected, data, L_mag, h_min_line, h_max_line, done_button, h_var_line)
    % Get current values
    min_val = get(src, 'Value');
    max_val = max(L_mag(:));
    
    % Update display
    caxis(h_plot.Parent, [min_val, max_val]);
    set(min_edit, 'String', sprintf('%.3f', min_val));
    set(h_min_line, 'Value', min_val);
    set(h_max_line, 'Value', max_val);
    
    % Find streaks
    streak_mask = L_mag >= min_val & L_mag <= max_val;
    [streak_rows, streak_cols] = find(streak_mask);
    
    % Correct streaks
    corrected = data;
    unique_cols = unique(streak_cols);
    
    for col = unique_cols'
        rows_in_col = streak_rows(streak_cols == col);
        
        % Process all streak points in this column
        if ~isempty(rows_in_col)
            % Find left neighbor
            left_col = col - 1;
            if left_col >= 1
                % Get left neighbor values
                left_vals = mean(data(rows_in_col, left_col));
                streak_avg = mean(data(rows_in_col, col));
                
                % Use left neighbor value
                corrected(rows_in_col, col) = data(rows_in_col, col) - streak_avg + left_vals;
            end
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
end


