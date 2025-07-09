function [corrected_data, streak_mask] = removeLocalStreaks(Y, slice_idx, min_value, max_streak_width, mode, auto)
%REMOVELOCALSTREAKS Interactive or automatic streak removal tool using multi-scale Laplacian detection
%   [corrected_data, streak_mask] = removeLocalStreaks(Y, slice_idx, min_value, max_streak_width, mode, auto)
%   If auto==true and min_value is provided, runs in non-interactive, non-visual mode.
%   Otherwise, runs interactively as before.
%
%   Inputs:
%       Y - 3D data array
%       slice_idx - Index of the slice to process (default: 150)
%       min_value - Optional minimum value for streak detection (default: interactive)
%       max_streak_width - Maximum streak width to detect (default: 3)
%       mode - 'valley', 'plateau', or 'both' (default: 'both')
%       auto - Boolean indicating whether to run in automatic mode (default: false)
%
%   Outputs:
%       corrected_data - The corrected image data
%       streak_mask - Binary mask indicating detected streaks
%
%   Example:
%       [corrected, mask] = removeLocalStreaks(Y, 150);
%       [corrected, mask] = removeLocalStreaks(Y, 150, 0.5, 5, 'valley');

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
if nargin < 5
    mode = 'both';
end
if nargin < 6
    auto = false;
end
corrected_data = [];
streak_mask = [];

% Process data
data = Y(:,:,slice_idx);
[rows, cols] = size(data);

% Compute global mean for fallback
global_mean = mean(data(:));

% Compute multi-scale Laplacians
L_valley = zeros(size(data));
L_plateau = zeros(size(data));
L_mag_combined = zeros(size(data));

for n = 1:max_streak_width
    % Shift left and right by n pixels
    data_left_n = [zeros(rows, n), data(:, 1:end-n)];
    data_right_n = [data(:, n+1:end), zeros(rows, n)];
    
    % Compute Laplacian at distance n
    L_n = data_left_n + data_right_n - 3*data;
    
    % Combine with previous scales (take maximum magnitude)
    L_valley = max(L_valley, L_n); % for valley (positive Laplacian)
    L_plateau = min(L_plateau, L_n); % for plateau (negative Laplacian)
    L_mag_n = abs(L_n);
    L_mag_combined = max(L_mag_combined, L_mag_n);
end

% Select Laplacian for mask and threshold range
if strcmpi(mode, 'valley')
    L_for_mask = L_valley;
    min_thr = 0;
    max_thr = max(L_for_mask(:));
    if min_thr == max_thr
        max_thr = min_thr + eps;
    end
    slider_range = [min_thr, max_thr];
elseif strcmpi(mode, 'plateau')
    L_for_mask = L_plateau;
    min_thr = min(L_for_mask(:));
    max_thr = 0;
    if min_thr == max_thr
        max_thr = min_thr + eps;
    elseif min_thr > max_thr
        tmp = min_thr; min_thr = max_thr; max_thr = tmp;
    end
    slider_range = [min_thr, max_thr];
else
    L_for_mask = L_mag_combined;
    min_thr = min(L_for_mask(:));
    max_thr = prctile(L_for_mask(:), 99.5);
    if min_thr == max_thr
        max_thr = min_thr + eps;
    end
    slider_range = [min_thr, max_thr];
end

% Set threshold value
if isempty(min_value)
    min_val = slider_range(1);
else
    min_val = max(slider_range(1), min(slider_range(2), min_value));
end

% Compute streak mask
if strcmpi(mode, 'valley')
    streak_mask = (L_for_mask > min_val);
elseif strcmpi(mode, 'plateau')
    streak_mask = (L_for_mask < min_val);
else
    streak_mask = (L_for_mask >= min_val);
end

% Compute correction
corrected = data;
[streak_rows, streak_cols] = find(streak_mask);
unique_cols = unique(streak_cols);
if strcmpi(mode, 'valley') || strcmpi(mode, 'plateau')
    if any(~streak_mask(:))
        unstreaked_mean = mean(data(~streak_mask));
    else
        unstreaked_mean = global_mean;
    end
    for col = unique_cols'
        rows_in_col = streak_rows(streak_cols == col);
        if ~isempty(rows_in_col)
            streak_avg = mean(data(rows_in_col, col));
            corrected(rows_in_col, col) = data(rows_in_col, col) - streak_avg + unstreaked_mean;
        end
    end
else
    for col = unique_cols'
        rows_in_col = streak_rows(streak_cols == col);
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
            if ~isempty(left_vals) && ~isempty(right_vals)
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
end
corrected_data = corrected;

% Visualization/UI only if not auto
if ~auto
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
    if strcmpi(mode, 'valley')
        hist_title = 'Histogram of Laplacian (valley mode)';
    elseif strcmpi(mode, 'plateau')
        hist_title = 'Histogram of Laplacian (plateau mode)';
    else
        hist_title = 'Histogram of Combined Laplacian Magnitude';
    end
    title(hist_title);
    axis square;
    hold on;
    hist_min = min(L_mag_combined(:));
    hist_max = max(L_mag_combined(:));
    h_min_line = xline(hist_min, 'r-', 'LineWidth', 2);
    h_max_line = xline(hist_max, 'r-', 'LineWidth', 2);
    ylabel('Count');
    xlim([hist_min, hist_max]);
    yyaxis right
    h_var_line = plot([hist_min, hist_max], [0, 0], 'b-', 'LineWidth', 2);
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
    initial_min = slider_range(1);
    if ~auto
        initial_min = min_value;
    end
    % Clamp initial_min to slider_range
    initial_min = max(slider_range(1), min(slider_range(2), initial_min));

    min_slider = uicontrol(panel, 'Style', 'slider', ...
        'Min', slider_range(1), ...
        'Max', slider_range(2), ...
        'Value', initial_min, ...
        'Position', [60, 5, 200, 20]);

    min_edit = uicontrol(panel, 'Style', 'edit', ...
        'String', sprintf('%.3f', initial_min), ...
        'Position', [270, 5, 80, 20]);

    done_button = uicontrol(panel, 'Style', 'pushbutton', ...
        'String', 'Done', ...
        'Position', [450, 5, 100, 40]);

    % Set up callbacks after all controls are created
    set(min_slider, 'Callback', @(src,event) updateContrast(src, event, h_mask, min_edit, h_corrected, data, L_for_mask, h_min_line, h_max_line, done_button, h_var_line, global_mean, mode));
    set(min_edit, 'Callback', @(src,event) updateFromText(src, event, min_slider, h_mask, h_corrected, data, L_for_mask, h_min_line, h_max_line, done_button, h_var_line, global_mean, mode));
    set(done_button, 'Callback', @(src,event) finish(src, h_corrected, h_mask));

    % Initialize display
    updateContrast(min_slider, [], h_mask, min_edit, h_corrected, data, L_for_mask, h_min_line, h_max_line, done_button, h_var_line, global_mean, mode);

    % Wait for the figure to be closed
    waitfor(h_fig);

    % Get results from the base workspace
    if evalin('base', 'exist(''temp_corrected_data'', ''var'')')
        corrected_data = evalin('base', 'temp_corrected_data');
        streak_mask = evalin('base', 'temp_streak_mask');
        evalin('base', 'clear temp_corrected_data temp_streak_mask');
    end
end

end

function updateFromText(src, ~, min_slider, h_mask, h_corrected, data, L_for_mask, h_min_line, h_max_line, done_button, h_var_line, global_mean, mode)
    % Get the value from the text input
    new_val = str2double(get(src, 'String'));
    
    % Validate the input
    if isnan(new_val)
        % If invalid input, reset to current slider value
        set(src, 'String', sprintf('%.3f', get(min_slider, 'Value')));
        return;
    end
    
    % Clamp the value to the valid range
    new_val = max(min(L_for_mask(:)), min(max(L_for_mask(:)), new_val));
    
    % Update the slider
    set(min_slider, 'Value', new_val);
    
    % Update the display
    updateContrast(min_slider, [], h_mask, src, h_corrected, data, L_for_mask, h_min_line, h_max_line, done_button, h_var_line, global_mean, mode);
    
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

function updateContrast(src, ~, h_mask, min_edit, h_corrected, data, L_for_mask, h_min_line, h_max_line, done_button, h_var_line, global_mean, mode)
    min_val = get(src, 'Value');
    max_val = get(src, 'Max');
    if min_val == max_val
        max_val = min_val + eps;
    elseif min_val > max_val
        tmp = min_val; min_val = max_val; max_val = tmp;
    end
    caxis(h_mask.Parent, [min_val, max_val]);
    set(min_edit, 'String', sprintf('%.3f', min_val));
    set(h_min_line, 'Value', min_val);
    set(h_max_line, 'Value', max_val);
    ax = ancestor(h_min_line, 'axes');
    hist_min = min(L_for_mask(:));
    hist_max = max(L_for_mask(:));
    if hist_min == hist_max
        hist_max = hist_min + eps;
    end
    xlim(ax, [hist_min, hist_max]);
    if strcmpi(mode, 'valley')
        streak_mask = (L_for_mask > min_val);
    elseif strcmpi(mode, 'plateau')
        streak_mask = (L_for_mask < min_val);
    else
        streak_mask = (L_for_mask >= min_val);
    end
    [streak_rows, streak_cols] = find(streak_mask);
    corrected = data;
    unique_cols = unique(streak_cols);
    if strcmpi(mode, 'valley') || strcmpi(mode, 'plateau')
        % Compute mean of unstreaked area
        if any(~streak_mask(:))
            unstreaked_mean = mean(data(~streak_mask));
        else
            unstreaked_mean = global_mean;
        end
        for col = unique_cols'
            rows_in_col = streak_rows(streak_cols == col);
            if ~isempty(rows_in_col)
                streak_avg = mean(data(rows_in_col, col));
                corrected(rows_in_col, col) = data(rows_in_col, col) - streak_avg + unstreaked_mean;
            end
        end
    else
        for col = unique_cols'
            rows_in_col = streak_rows(streak_cols == col);
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
                if ~isempty(left_vals) && ~isempty(right_vals)
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
    end
    var_corrected = var(corrected(:));
    set(h_var_line, 'YData', [var_corrected, var_corrected]);
    set(h_corrected, 'CData', corrected);
    caxis(h_corrected.Parent, [min(data(:)), max(data(:))]);
    set(done_button, 'UserData', struct('min_val', min_val, 'max_val', max_val));
    drawnow;
    set(h_mask, 'CData', streak_mask);
end


