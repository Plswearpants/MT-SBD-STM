function [Y_masked, defect_mask] = thresholdDefects(Y, threshold_slice)
    % Interactive thresholding of defects in 3D data
    %
    % Inputs:
    %   Y: 3D data array [height x width x depth]
    %   threshold_slice: (optional) slice number to use for thresholding
    %                    if not provided, will ask user to select
    %
    % Outputs:
    %   Y_masked: 3D data with defects masked/interpolated
    %   defect_mask: binary mask of defects (true for defects)
    
    % Create a figure for thresholding
    f_threshold = figure('Name', 'Defect Thresholding', 'Position', [100, 100, 1200, 600]);
    
    % Create panels for image, histogram, and controls
    img_panel = uipanel('Parent', f_threshold, 'Position', [0.05, 0.1, 0.4, 0.8]);
    hist_panel = uipanel('Parent', f_threshold, 'Position', [0.5, 0.1, 0.45, 0.8]);
    
    % If threshold_slice not provided, ask user to select
    if nargin < 2 || isempty(threshold_slice)
        f_select = figure;
        d3gridDisplay(Y, 'dynamic');
        threshold_slice = input('Enter slice number for thresholding: ');
        close(f_select);
    end
    
    % Get the selected slice data
    slice_data = Y(:,:,threshold_slice);
    
    % Calculate initial min and max values
    data_min = min(slice_data(:));
    data_max = max(slice_data(:));
    
    % Initialize threshold values
    lower_threshold = data_min;
    upper_threshold = data_max;
    
    % Create axes for image and histogram
    img_ax = axes('Parent', img_panel);
    hist_ax = axes('Parent', hist_panel);
    
    % Function to update the display
    function updateDisplay()
        % Update image with thresholding
        cla(img_ax);
        axes(img_ax);
        masked_data = slice_data;
        mask = (slice_data < lower_threshold) | (slice_data > upper_threshold);
        imagesc(masked_data);
        colormap(img_ax, 'jet');
        colorbar(img_ax);
        title(img_ax, sprintf('Slice %d with Threshold [%.2f, %.2f]', threshold_slice, lower_threshold, upper_threshold));
        axis(img_ax, 'image');
        
        % Apply autocontrast
        caxis(img_ax, [min(masked_data(:)), max(masked_data(:))]);
        
        % Highlight masked regions
        hold(img_ax, 'on');
        h = imagesc(img_ax, ones(size(masked_data)));
        set(h, 'AlphaData', 0.5*double(mask));
        colormap(img_ax, 'jet');
        hold(img_ax, 'off');
        
        % Update histogram
        cla(hist_ax);
        axes(hist_ax);
        histogram(hist_ax, slice_data(:), 100);
        hold(hist_ax, 'on');
        yl = ylim(hist_ax);
        plot(hist_ax, [lower_threshold, lower_threshold], [0, yl(2)], 'r-', 'LineWidth', 2);
        plot(hist_ax, [upper_threshold, upper_threshold], [0, yl(2)], 'r-', 'LineWidth', 2);
        title(hist_ax, 'Histogram with Threshold Lines');
        xlabel(hist_ax, 'Pixel Value');
        ylabel(hist_ax, 'Frequency');
        hold(hist_ax, 'off');
    end
    
    % Create sliders for thresholding
    slider_lower = uicontrol('Parent', f_threshold, 'Style', 'slider', ...
        'Position', [100, 40, 400, 20], ...
        'Min', data_min, 'Max', data_max, 'Value', data_min, ...
        'Callback', @(src, ~) updateLowerThreshold(src.Value));
    
    slider_upper = uicontrol('Parent', f_threshold, 'Style', 'slider', ...
        'Position', [600, 40, 400, 20], ...
        'Min', data_min, 'Max', data_max, 'Value', data_max, ...
        'Callback', @(src, ~) updateUpperThreshold(src.Value));
    
    % Create text labels for sliders
    uicontrol('Parent', f_threshold, 'Style', 'text', ...
        'Position', [100, 60, 100, 20], ...
        'String', 'Lower Threshold:');
    
    uicontrol('Parent', f_threshold, 'Style', 'text', ...
        'Position', [600, 60, 100, 20], ...
        'String', 'Upper Threshold:');
    
    % Create text displays for current threshold values
    lower_text = uicontrol('Parent', f_threshold, 'Style', 'text', ...
        'Position', [200, 60, 100, 20], ...
        'String', num2str(lower_threshold));
    
    upper_text = uicontrol('Parent', f_threshold, 'Style', 'text', ...
        'Position', [700, 60, 100, 20], ...
        'String', num2str(upper_threshold));
    
    % Create confirm button
    confirm_btn = uicontrol('Parent', f_threshold, 'Style', 'pushbutton', ...
        'Position', [500, 10, 100, 30], ...
        'String', 'Confirm', ...
        'Callback', @confirmThreshold);
    
    % Callback functions
    function updateLowerThreshold(value)
        lower_threshold = value;
        if lower_threshold > upper_threshold
            lower_threshold = upper_threshold;
            slider_lower.Value = lower_threshold;
        end
        lower_text.String = num2str(lower_threshold);
        updateDisplay();
    end
    
    function updateUpperThreshold(value)
        upper_threshold = value;
        if upper_threshold < lower_threshold
            upper_threshold = lower_threshold;
            slider_upper.Value = upper_threshold;
        end
        upper_text.String = num2str(upper_threshold);
        updateDisplay();
    end
    
    function confirmThreshold(~, ~)
        % Apply thresholding to all slices
        fprintf('Applying threshold [%.2f, %.2f] to all slices...\n', lower_threshold, upper_threshold);
        
        % Create a mask for defects across all slices
        defect_mask = false(size(Y));
        
        % Apply thresholding to each slice
        for i = 1:size(Y, 3)
            defect_mask(:,:,i) = (Y(:,:,i) < lower_threshold) | (Y(:,:,i) > upper_threshold);
        end
        
        % Initialize masked data
        Y_masked = Y;
        
        % Replace values outside threshold with limit values
        Y_masked(Y < lower_threshold) = lower_threshold;
        Y_masked(Y > upper_threshold) = upper_threshold;
        
        % Display the result
        figure;
        subplot(1,2,1);
        imagesc(Y(:,:,threshold_slice));
        colorbar;
        title('Original Data');
        axis image;
        caxis([min(Y(:)), max(Y(:))]);  % Apply autocontrast
        
        subplot(1,2,2);
        imagesc(Y_masked(:,:,threshold_slice));
        colorbar;
        title('After Defect Removal');
        axis image;
        caxis([min(Y_masked(:)), max(Y_masked(:))]);  % Apply autocontrast
        
        % Close the thresholding figure
        close(f_threshold);
    end
    
    % Initial display
    updateDisplay();
    
    % Wait for user to confirm
    uiwait(f_threshold);
end 