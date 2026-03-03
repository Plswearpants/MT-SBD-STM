function [corrected_data, streak_mask, final_min_val] = removeLocalStreaks(Y, slice_idx, min_value, max_streak_width, mode, auto)
%REMOVELOCALSTREAKS Interactive or automatic streak removal tool using multi-scale Laplacian detection
%   [corrected_data, streak_mask] = removeLocalStreaks(Y, slice_idx, min_value, max_streak_width, mode, auto)
%   [corrected_data, streak_mask, final_min_val] = removeLocalStreaks(...)
%   When auto==false (interactive), on Done final_min_val is the chosen threshold so the
%   caller can compute effective factor as final_min_val/min_low for reuse or recording.
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
final_min_val = [];

% Process data
data = Y(:,:,slice_idx);
[rows, cols] = size(data);

% Laplacian and threshold range (shared with streak_correction via streakCore)
[L_for_mask, slider_range] = streakCore('laplacian', data, mode, max_streak_width);
% Population-ranked threshold list for UI slider (linear in rank, not value)
lower_list_ui = streakCore('threshold_list', L_for_mask, 100);

if isempty(min_value)
    min_val = slider_range(1);
else
    min_val = max(slider_range(1), min(slider_range(2), min_value));
end

streak_mask = streakCore('mask', L_for_mask, mode, min_val);
if strcmpi(mode, 'valley') || strcmpi(mode, 'plateau')
    corrected_data = streakCore('correct', data, streak_mask, 'background_mean');
else
    corrected_data = streakCore('correct', data, streak_mask, 'neighbor_interp');
end
global_mean = mean(data(:));  % for UI callbacks

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

    % Plot binary mask (0/1: detected streaks)
    subplot(2, 2, 2);
    h_mask = imagesc(zeros(size(L_for_mask)));
    title('Detected Streak Mask');
    axis square;
    colormap(gca, gray);
    caxis(h_mask.Parent, [0, 1]);  % binary mask: 0 = no streak, 1 = streak
    colorbar;

    % Plot histogram of the same Laplacian used for the mask (so threshold line matches)
    subplot(2, 2, 3);
    yyaxis left
    h_hist = histogram(L_for_mask);
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
    hist_min = min(L_for_mask(:));
    hist_max = max(L_for_mask(:));
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

    % Initialize slider position from input (maps to nearest population-ranked threshold)
    initial_min = slider_range(1);
    if ~auto && ~isempty(min_value)
        initial_min = min_value;
    end
    % Clamp initial_min to slider_range
    initial_min = max(slider_range(1), min(slider_range(2), initial_min));

    % Slider operates over indices into lower_list_ui (linear in population rank)
    n_levels = numel(lower_list_ui);
    if n_levels < 1
        lower_list_ui = slider_range(1);
        n_levels = 1;
    end
    [~, initial_idx] = min(abs(lower_list_ui - initial_min));
    min_slider = uicontrol(panel, 'Style', 'slider', ...
        'Min', 1, ...
        'Max', n_levels, ...
        'Value', initial_idx, ...
        'SliderStep', [1/max(1,n_levels-1), min(1,10/max(1,n_levels-1))], ...
        'Position', [60, 5, 200, 20]);

    min_edit = uicontrol(panel, 'Style', 'edit', ...
        'String', sprintf('%.3f', lower_list_ui(initial_idx)), ...
        'Position', [270, 5, 80, 20]);

    done_button = uicontrol(panel, 'Style', 'pushbutton', ...
        'String', 'Done', ...
        'Position', [450, 5, 100, 40]);

    % Set up callbacks after all controls are created
    set(min_slider, 'Callback', @(src,event) updateContrast(src, event, h_mask, min_edit, h_corrected, data, L_for_mask, lower_list_ui, h_min_line, h_max_line, done_button, h_var_line, global_mean, mode));
    set(min_edit, 'Callback', @(src,event) updateFromText(src, event, min_slider, h_mask, h_corrected, data, L_for_mask, lower_list_ui, h_min_line, h_max_line, done_button, h_var_line, global_mean, mode));
    set(done_button, 'Callback', @(src,event) finish(src, h_corrected, h_mask));

    % Initialize display
    updateContrast(min_slider, [], h_mask, min_edit, h_corrected, data, L_for_mask, lower_list_ui, h_min_line, h_max_line, done_button, h_var_line, global_mean, mode);

    % Wait for the figure to be closed
    waitfor(h_fig);

    % Get results from the base workspace (including final slider threshold for factor feedback)
    if evalin('base', 'exist(''temp_corrected_data'', ''var'')')
        corrected_data = evalin('base', 'temp_corrected_data');
        streak_mask = evalin('base', 'temp_streak_mask');
        if evalin('base', 'exist(''temp_final_min_val_remove'', ''var'')')
            final_min_val = evalin('base', 'temp_final_min_val_remove');
            evalin('base', 'clear temp_final_min_val_remove');
        end
        evalin('base', 'clear temp_corrected_data temp_streak_mask');
    end
end

end

function updateFromText(src, ~, min_slider, h_mask, h_corrected, data, L_for_mask, lower_list_ui, h_min_line, h_max_line, done_button, h_var_line, global_mean, mode)
    % Get the value from the text input
    new_val = str2double(get(src, 'String'));
    
    % Validate the input
    if isnan(new_val)
        % If invalid input, reset to current slider-based threshold
        idx = round(get(min_slider, 'Value'));
        idx = max(1, min(numel(lower_list_ui), idx));
        set(src, 'String', sprintf('%.3f', lower_list_ui(idx)));
        return;
    end
    
    % Map typed threshold to nearest population-ranked level
    if isempty(lower_list_ui)
        return;
    end
    [~, idx] = min(abs(lower_list_ui - new_val));
    set(min_slider, 'Value', idx);
    
    % Update the display (this will also refresh the text with the snapped value)
    updateContrast(min_slider, [], h_mask, src, h_corrected, data, L_for_mask, lower_list_ui, h_min_line, h_max_line, done_button, h_var_line, global_mean, mode);
    
    % Force immediate update of all figures
    drawnow;
end

function finish(src, h_corrected, h_mask)
    % Get the current values from the button's UserData
    user_data = get(src, 'UserData');
    min_val = user_data(1).min_val;
    max_val = user_data(1).max_val;
    
    % Store final slider threshold for caller (effective factor = min_val / min_low)
    assignin('base', 'temp_final_min_val_remove', min_val);
    
    % CData of h_mask is the binary streak mask (0/1) we set in updateContrast
    corrected_data = get(h_corrected, 'CData');
    streak_mask = logical(get(h_mask, 'CData'));
    
    % Store results in base workspace
    assignin('base', 'temp_corrected_data', corrected_data);
    assignin('base', 'temp_streak_mask', streak_mask);
    
    % Close the figure
    close(gcf);
end

function updateContrast(src, ~, h_mask, min_edit, h_corrected, data, L_for_mask, lower_list_ui, h_min_line, h_max_line, done_button, h_var_line, ~, mode)
    % Slider value is an index into lower_list_ui (population-ranked thresholds)
    idx = round(get(src, 'Value'));
    idx = max(1, min(numel(lower_list_ui), idx));
    set(src, 'Value', idx);
    min_val = lower_list_ui(idx);
    set(min_edit, 'String', sprintf('%.3f', min_val));
    set(h_min_line, 'Value', min_val);
    set(h_max_line, 'Value', min_val);
    ax = ancestor(h_min_line, 'axes');
    hist_min = min(L_for_mask(:));
    hist_max = max(L_for_mask(:));
    if hist_min == hist_max
        hist_max = hist_min + eps;
    end
    xlim(ax, [hist_min, hist_max]);
    streak_mask = streakCore('mask', L_for_mask, mode, min_val);
    if strcmpi(mode, 'valley') || strcmpi(mode, 'plateau')
        corrected = streakCore('correct', data, streak_mask, 'background_mean');
    else
        corrected = streakCore('correct', data, streak_mask, 'neighbor_interp');
    end
    set(h_var_line, 'YData', [var(corrected(:)), var(corrected(:))]);
    set(h_corrected, 'CData', corrected);
    caxis(h_corrected.Parent, [min(data(:)), max(data(:))]);
    set(done_button, 'UserData', struct('min_val', min_val, 'max_val', max(L_for_mask(:))));
    set(h_mask, 'CData', double(streak_mask));
    caxis(h_mask.Parent, [0, 1]);
    drawnow;
end


