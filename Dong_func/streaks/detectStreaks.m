function detectStreaks(I, threshold, min_length)
% detectStreaks - Detect horizontal streaks (vertical edges) in STM image
% I          : Input 2D STM image
% threshold  : Threshold for detecting streaks (e.g., 10)
% min_length : Minimum length (in pixels) of streaks to retain

if nargin < 2
    threshold = 10;
end
if nargin < 3
    min_length = 10;
end

% Convert to double for filtering
I = double(I);

% Vertical gradient kernel (for horizontal streaks)
K = [-1; 0; 1];

% Apply 2D correlation
edgeResponse = conv2(I, K, 'same');

% Create figure for interactive visualization
h_fig = figure('Name', 'Interactive Streak Detection', 'Position', [100, 100, 1200, 800]);

% Plot original image
subplot(2, 2, 1);
h_orig = imagesc(I);
title('Original Image');
axis square;
colormap gray;
colorbar;

% Plot edge response
subplot(2, 2, 2);
h_edge = imagesc(abs(edgeResponse));
title('Edge Response Magnitude');
axis square;
colormap gray;
colorbar;

% Plot histogram
subplot(2, 2, 3);
h_hist = histogram(abs(edgeResponse(:)));
title('Histogram of Edge Response');
xlabel('Edge Response Magnitude');
ylabel('Count');
hold on;
h_thresh_line = xline(threshold, 'r-', 'LineWidth', 2, 'Label', 'Threshold');
hold off;

% Plot detected streaks
subplot(2, 2, 4);
h_streak = imagesc(I);
title('Detected Streaks');
axis square;
colormap gray;
colorbar;
hold on;

% Create controls panel
panel = uipanel('Position', [0.1, 0.02, 0.8, 0.05]);

% Threshold slider
thresh_label = uicontrol(panel, 'Style', 'text', 'String', 'Threshold:', 'Position', [10, 5, 60, 20]);
thresh_slider = uicontrol(panel, 'Style', 'slider', ...
    'Min', min(abs(edgeResponse(:))), ...
    'Max', max(abs(edgeResponse(:))), ...
    'Value', threshold, ...
    'Position', [80, 5, 200, 20]);
thresh_edit = uicontrol(panel, 'Style', 'edit', ...
    'String', sprintf('%.1f', threshold), ...
    'Position', [290, 5, 80, 20]);

% Min length slider
length_label = uicontrol(panel, 'Style', 'text', 'String', 'Min Length:', 'Position', [380, 5, 70, 20]);
length_slider = uicontrol(panel, 'Style', 'slider', ...
    'Min', 1, ...
    'Max', 50, ...
    'Value', min_length, ...
    'Position', [460, 5, 150, 20]);
length_edit = uicontrol(panel, 'Style', 'edit', ...
    'String', sprintf('%d', min_length), ...
    'Position', [620, 5, 60, 20]);

% Done button
done_button = uicontrol(panel, 'Style', 'pushbutton', ...
    'String', 'Done', ...
    'Position', [690, 5, 100, 40]);

% Set up callbacks
set(thresh_slider, 'Callback', @(src,event) updateDetection(src, event, thresh_edit, length_slider, length_edit, h_thresh_line, h_streak, I, edgeResponse));
set(thresh_edit, 'Callback', @(src,event) updateFromText(src, event, thresh_slider, length_slider, length_edit, h_thresh_line, h_streak, I, edgeResponse));
set(length_slider, 'Callback', @(src,event) updateDetection(thresh_slider, event, thresh_edit, src, length_edit, h_thresh_line, h_streak, I, edgeResponse));
set(length_edit, 'Callback', @(src,event) updateFromTextLength(src, event, thresh_slider, thresh_edit, length_slider, h_thresh_line, h_streak, I, edgeResponse));
set(done_button, 'Callback', @(src,event) finish(src, h_streak));

% Initialize display
updateDetection(thresh_slider, [], thresh_edit, length_slider, length_edit, h_thresh_line, h_streak, I, edgeResponse);

% Wait for the figure to be closed
waitfor(h_fig);

end

function updateFromText(src, ~, thresh_slider, length_slider, length_edit, h_thresh_line, h_streak, I, edgeResponse)
    % Get the value from the text input
    new_val = str2double(get(src, 'String'));
    
    % Validate the input
    if isnan(new_val)
        % If invalid input, reset to current slider value
        set(src, 'String', sprintf('%.1f', get(thresh_slider, 'Value')));
        return;
    end
    
    % Clamp the value to the valid range
    new_val = max(min(abs(edgeResponse(:))), min(max(abs(edgeResponse(:))), new_val));
    
    % Update the slider
    set(thresh_slider, 'Value', new_val);
    
    % Update the display
    updateDetection(thresh_slider, [], src, length_slider, length_edit, h_thresh_line, h_streak, I, edgeResponse);
    
    % Force immediate update
    drawnow;
    
    % Update the text field with the clamped value if it was changed
    if new_val ~= str2double(get(src, 'String'))
        set(src, 'String', sprintf('%.1f', new_val));
    end
end

function updateFromTextLength(src, ~, thresh_slider, thresh_edit, length_slider, h_thresh_line, h_streak, I, edgeResponse)
    % Get the value from the text input
    new_val = str2double(get(src, 'String'));
    
    % Validate the input
    if isnan(new_val) || new_val < 1
        % If invalid input, reset to current slider value
        set(src, 'String', sprintf('%d', get(length_slider, 'Value')));
        return;
    end
    
    % Clamp the value to the valid range
    new_val = max(1, min(50, new_val));
    
    % Update the slider
    set(length_slider, 'Value', new_val);
    
    % Update the display
    updateDetection(thresh_slider, [], thresh_edit, length_slider, src, h_thresh_line, h_streak, I, edgeResponse);
    
    % Force immediate update
    drawnow;
    
    % Update the text field with the clamped value if it was changed
    if new_val ~= str2double(get(src, 'String'))
        set(src, 'String', sprintf('%d', new_val));
    end
end

function finish(src, h_streak)
    % Get the current streak map from the plot
    streak_map = get(h_streak, 'UserData');
    
    % Store results in base workspace
    assignin('base', 'detected_streaks', streak_map);
    
    % Close the figure
    close(gcf);
end

function updateDetection(thresh_src, ~, thresh_edit, length_src, length_edit, h_thresh_line, h_streak, I, edgeResponse)
    % Get current values
    threshold = get(thresh_src, 'Value');
    min_length = get(length_src, 'Value');
    
    % Update threshold line
    set(h_thresh_line, 'Value', threshold);
    set(h_thresh_line, 'Label', sprintf('Threshold: %.1f', threshold));
    
    % Update text fields
    set(thresh_edit, 'String', sprintf('%.1f', threshold));
    set(length_edit, 'String', sprintf('%d', min_length));
    
    % Binary thresholding
    streakMap = abs(edgeResponse) > threshold;
    
    % Remove small objects (denoise)
    streakMap = bwareaopen(streakMap, min_length);
    
    % Update streak visualization
    set(h_streak, 'CData', I);
    hold(h_streak.Parent, 'on');
    
    % Clear previous boundaries
    children = get(h_streak.Parent, 'Children');
    for i = 1:length(children)
        if strcmp(get(children(i), 'Type'), 'line')
            delete(children(i));
        end
    end
    
    % Draw new boundaries
    [B,~] = bwboundaries(streakMap, 'noholes');
    for k = 1:length(B)
        boundary = B{k};
        plot(h_streak.Parent, boundary(:,2), boundary(:,1), 'r', 'LineWidth', 1);
    end
    hold(h_streak.Parent, 'off');
    
    % Store streak map in the plot for later retrieval
    set(h_streak, 'UserData', streakMap);
    
    % Force immediate update
    drawnow;
end