function [corrected_data, streak_mask, streak_indices, final_min_val] = interpolateLocalStreaks(Y, slice_idx, min_value, provided_streak_indices, auto)
%INTERPOLATELOCALSTREAKS Interactive or automatic streak interpolating tool using Laplacian and neighbor interpolation
%   Uses streakCore for Laplacian (mode 'both', max_streak_width 1), mask, and horizontal_avg correction.
%   [corrected_data, streak_mask, streak_indices] = interpolateLocalStreaks(...)
%   [corrected_data, streak_mask, streak_indices, final_min_val] = interpolateLocalStreaks(...)
%   When interactive, on Done final_min_val is the chosen threshold for factor feedback (effective factor = final_min_val/min_low).
%   If auto==true and min_value is provided, runs in non-interactive, non-visual mode.
%   Otherwise, runs interactively as before.

if nargin < 2
    slice_idx = 150;
end
if nargin < 3
    min_value = [];
end
if nargin < 4
    provided_streak_indices = [];
end
if nargin < 5
    auto = false;
end
corrected_data = [];
streak_mask = [];
streak_indices = [];
final_min_val = [];

data = single(Y(:,:,slice_idx));
[rows, cols] = size(data);

% If streak indices are provided, apply horizontal_avg on that mask (same logic as core)
if nargin >= 4 && ~isempty(provided_streak_indices)
    streak_mask = false(size(data));
    streak_mask(sub2ind(size(data), provided_streak_indices(:,1), provided_streak_indices(:,2))) = true;
    corrected_data = streakCore('correct', double(data), streak_mask, 'horizontal_avg');
    valid_streaks = streak_mask;
    valid_streaks(:,1) = false;
    valid_streaks(:,end) = false;
    valid_streaks = valid_streaks & circshift(streak_mask, [0 -1]) & circshift(streak_mask, [0 1]);
    streak_mask = valid_streaks;
    [streak_rows, streak_cols] = find(streak_mask);
    streak_indices = [streak_rows, streak_cols];
    return;
end

% Single-scale Laplacian (n=1) and horizontal_avg via streakCore
[L_for_mask, slider_range] = streakCore('laplacian', double(data), 'both', 1);
L_mag = L_for_mask;
% Population-ranked threshold list for UI slider (linear in rank, not value)
lower_list_ui = streakCore('threshold_list', L_mag, 100);
if isempty(min_value)
    min_val = slider_range(1);
else
    min_val = min_value;
end
streak_mask_raw = streakCore('mask', L_for_mask, 'both', min_val);
corrected_data = streakCore('correct', double(data), streak_mask_raw, 'horizontal_avg');
% Pixels actually interpolated = middle cols with both neighbors in raw mask
valid_streaks = streak_mask_raw;
valid_streaks(:,1) = false;
valid_streaks(:,end) = false;
valid_streaks = valid_streaks & circshift(streak_mask_raw, [0 -1]) & circshift(streak_mask_raw, [0 1]);
streak_mask = valid_streaks;
[streak_rows, streak_cols] = find(streak_mask);
streak_indices = [streak_rows, streak_cols];

% Visualization/UI only if not auto and no provided_streak_indices
if ~auto && isempty(provided_streak_indices)
    % --- Interactive UI as before ---
    h_fig = figure('Name', 'X-Direction Laplacian Streak Removal Analysis', 'Position', [100, 100, 1200, 800]);
    subplot(2, 2, 1);
    h_orig = imagesc(data);
    title('Original Data');
    axis square;
    colormap gray;
    colorbar;
    subplot(2, 2, 2);
    h_mask_plot = imagesc(false(size(L_mag)));
    title('Detected Streak Mask (used for interpolation)');
    axis square;
    colormap(gca, gray);
    caxis(gca, [0, 1]);
    colorbar;
    subplot(2, 2, 3);
    h_hist = histogram(L_mag);
    title('Histogram of Laplacian Magnitude');
    axis square;
    hold on;
    h_min_line = xline(min(L_mag(:)), 'r-', 'LineWidth', 2);
    h_max_line = xline(max(L_mag(:)), 'r-', 'LineWidth', 2);
    xlim([min(L_mag(:)), max(L_mag(:))/2]);
    hold off;
    subplot(2, 2, 4);
    h_corrected = imagesc(data);
    title('Corrected Image');
    axis square;
    colormap gray;
    colorbar;
    panel = uipanel('Position', [0.1, 0.02, 0.8, 0.05]);
    initial_min = min(L_mag(:));
    if ~isempty(min_value)
        initial_min = min_value;
    end
    uicontrol(panel, 'Style', 'text', 'String', 'Min:', 'Position', [10, 5, 40, 20]);
    % Slider operates over indices into lower_list_ui (linear in population rank)
    n_levels = numel(lower_list_ui);
    if n_levels < 1
        lower_list_ui = initial_min;
        n_levels = 1;
    end
    [~, initial_idx] = min(abs(lower_list_ui - initial_min));
    min_slider = uicontrol(panel, 'Style', 'slider', ...
        'Min', 1, ...
        'Max', n_levels, ...
        'Value', initial_idx, ...
        'SliderStep', [1/max(1,n_levels-1), min(1,10/max(1,n_levels-1))], ...
        'Position', [60, 5, 300, 20]);
    min_text = uicontrol(panel, 'Style', 'text', ...
        'String', sprintf('%.3f', lower_list_ui(initial_idx)), ...
        'Position', [370, 5, 60, 20]);
    done_button = uicontrol(panel, 'Style', 'pushbutton', ...
        'String', 'Done', ...
        'Position', [450, 5, 100, 40], ...
        'Callback', @(src,event) finish(src, h_corrected, h_mask_plot));
    set(min_slider, 'Callback', @(src,event) debouncedUpdate(src, event, h_mask_plot, min_text, h_corrected, data, L_mag, lower_list_ui, h_min_line, h_max_line, done_button));
    updateContrast(min_slider, [], h_mask_plot, min_text, h_corrected, data, L_mag, lower_list_ui, h_min_line, h_max_line, done_button);
    waitfor(h_fig);
    if evalin('base', 'exist(''temp_corrected_data'', ''var'')')
        corrected_data = evalin('base', 'temp_corrected_data');
        streak_mask = evalin('base', 'temp_streak_mask');
        streak_indices = evalin('base', 'temp_streak_indices');
        if evalin('base', 'exist(''temp_final_min_val_interp'', ''var'')')
            final_min_val = evalin('base', 'temp_final_min_val_interp');
            evalin('base', 'clear temp_final_min_val_interp');
        end
        evalin('base', 'clear temp_corrected_data temp_streak_mask temp_streak_indices');
    end
end

end

function finish(src, h_corrected, h_mask)
    % Get the current values from the button's UserData
    user_data = get(src, 'UserData');
    min_val = user_data(1).min_val;
    max_val = user_data(1).max_val;
    
    % Store final slider threshold for caller (effective factor = min_val / min_low)
    assignin('base', 'temp_final_min_val_interp', min_val);
    
    % CData is the mask we displayed (valid_streaks) and used for correction
    corrected_data = get(h_corrected, 'CData');
    streak_mask = logical(get(h_mask, 'CData'));
    [streak_rows, streak_cols] = find(streak_mask);
    streak_indices = [streak_rows, streak_cols];
    
    % Store results in base workspace
    assignin('base', 'temp_corrected_data', corrected_data);
    assignin('base', 'temp_streak_mask', streak_mask);
    assignin('base', 'temp_streak_indices', streak_indices);
    
    % Close the figure
    close(gcf);
end

function debouncedUpdate(src, event, h_mask, min_text, h_corrected, data, L_mag, lower_list_ui, h_min_line, h_max_line, done_button)
    persistent lastUpdate
    if isempty(lastUpdate)
        lastUpdate = tic;
    end
    
    % Only update if 0.1 seconds have passed since last update
    if toc(lastUpdate) > 0.1
        updateContrast(src, event, h_mask, min_text, h_corrected, data, L_mag, lower_list_ui, h_min_line, h_max_line, done_button);
        lastUpdate = tic;
    end
end

function updateContrast(src, ~, h_mask, min_text, h_corrected, data, L_mag, lower_list_ui, h_min_line, h_max_line, done_button)
    % Slider value is an index into lower_list_ui (population-ranked thresholds)
    idx = round(get(src, 'Value'));
    idx = max(1, min(numel(lower_list_ui), idx));
    set(src, 'Value', idx);
    min_val = lower_list_ui(idx);

    streak_mask_raw = streakCore('mask', L_mag, 'both', min_val);
    corrected = streakCore('correct', double(data), streak_mask_raw, 'horizontal_avg');
    valid_streaks = streak_mask_raw;
    valid_streaks(:,1) = false;
    valid_streaks(:,end) = false;
    valid_streaks = valid_streaks & circshift(streak_mask_raw, [0 -1]) & circshift(streak_mask_raw, [0 1]);
    set(h_mask, 'CData', double(valid_streaks));
    caxis(h_mask.Parent, [0, 1]);
    caxis(h_corrected.Parent, [min(data(:)), max(data(:))]);
    set(min_text, 'String', sprintf('%.3f', min_val));
    set(h_min_line, 'Value', min_val);
    set(h_max_line, 'Value', max(L_mag(:)));
    set(h_corrected, 'CData', corrected);
    set(done_button, 'UserData', struct('min_val', min_val, 'max_val', max(L_mag(:))));
end


