function fit_polynomial_surface_ui(data)
% FIT_POLYNOMIAL_SURFACE_UI Interactive UI for polynomial surface fitting
%   fit_polynomial_surface_ui(data)
%   - data: 2D array to fit
%   Allows interactive mask drawing and live adjustment of polynomial degrees.

if ndims(data) ~= 2
    error('This UI only supports 2D data. For 3D, use fit_polynomial_surface.');
end

% Initial degrees
deg_x = 2;
deg_y = 2;

% Create UI figure
fig = uifigure('Name', 'Polynomial Surface Fitting UI', 'Position', [100 100 1400 600]);

% Axes for images
ax1 = uiaxes(fig, 'Position', [30 220 350 350]);
ax2 = uiaxes(fig, 'Position', [400 220 350 350]);
ax3 = uiaxes(fig, 'Position', [770 220 350 350]);

% Degree sliders
lblX = uilabel(fig, 'Position', [30 170 120 22], 'Text', 'Degree X:');
sliderX = uislider(fig, 'Position', [120 180 200 3], 'Limits', [0 10], 'Value', deg_x, 'MajorTicks', 0:10);
lblY = uilabel(fig, 'Position', [30 130 120 22], 'Text', 'Degree Y:');
sliderY = uislider(fig, 'Position', [120 140 200 3], 'Limits', [0 10], 'Value', deg_y, 'MajorTicks', 0:10);

% Mask drawing button
btnDraw = uibutton(fig, 'push', 'Text', 'Draw Mask(s)', 'Position', [400 140 120 40]);

% Info label
infoLbl = uilabel(fig, 'Position', [30 80 1000 22], 'Text', 'Draw mask(s) first, then adjust degrees.');

% Store mask data in app
appdata.masks = {};

% Draw initial fit
update_fit();

% Callbacks
sliderX.ValueChangedFcn = @(src, event) onDegreeChange();
sliderY.ValueChangedFcn = @(src, event) onDegreeChange();
btnDraw.ButtonPushedFcn = @(src, event) onDrawMask();

    function onDegreeChange()
        deg_x = round(sliderX.Value);
        deg_y = round(sliderY.Value);
        update_fit();
    end

    function onDrawMask()
        masks = draw_interactive_masks(data, ax1);
        appdata.masks = masks;
        btnDraw.Enable = 'off'; % Disable the button after use
        update_fit();
    end

    function update_fit()
        % Use current degrees and masks
        degs = [deg_x deg_y];
        if isempty(appdata.masks)
            masks_to_use = {};
        else
            masks_to_use = appdata.masks;
        end
        [fitted, subtracted] = fit_polynomial_surface(data, degs, 'Mask', masks_to_use);
        % Plot original with masks
        cla(ax1);
        imagesc(ax1, data);
        colormap(ax1, 'gray');
        axis(ax1, 'equal', 'tight');
        title(ax1, 'Original Data + Masks');
        hold(ax1, 'on');
        for i = 1:length(masks_to_use)
            plot(ax1, masks_to_use{i}(:,1), masks_to_use{i}(:,2), 'r-', 'LineWidth', 2);
        end
        hold(ax1, 'off');
        % Plot fitted
        cla(ax2);
        imagesc(ax2, fitted);
        colormap(ax2, 'gray');
        axis(ax2, 'equal', 'tight');
        title(ax2, sprintf('Fitted Surface (deg_x=%d, deg_y=%d)', deg_x, deg_y));
        % Plot subtracted
        cla(ax3);
        imagesc(ax3, subtracted);
        colormap(ax3, 'gray');
        axis(ax3, 'equal', 'tight');
        title(ax3, 'Subtracted Surface');
    end
end

function masks = draw_interactive_masks(data, ax)
    % Draw masks interactively on the given axes
    figureHandle = ancestor(ax, 'figure');
    figure(figureHandle); % bring to front
    axes(ax); % set current axes
    masks = {};
    while true
        h = drawpolygon('Color', 'r');
        if isempty(h.Position)
            break;
        end
        masks{end+1} = h.Position;
        h.Color = 'g';
        % Prompt for another?
        answer = questdlg('Add another mask?', 'Mask', 'Yes', 'No', 'No');
        if strcmp(answer, 'No')
            break;
        end
    end
end 