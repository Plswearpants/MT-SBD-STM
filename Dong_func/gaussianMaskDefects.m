function [Y_masked, defect_mask] = gaussianMaskDefects(Y, threshold_slice)
    % Apply Gaussian masking to defects in 3D data
    % g_M(r,E) = g(r,E) × (1 - M(r - r_0, σ))
    % M is a truncated Gaussian with maximum value = 0.99
    % σ is approximately half width of the defect center
    %
    % Inputs:
    %   Y: 3D data array [height x width x depth]
    %   threshold_slice: (optional) slice number to use for defect selection
    %                    if not provided, will ask user to select
    %
    % Outputs:
    %   Y_masked: 3D data with defects masked using Gaussian suppression
    %   defect_mask: binary mask showing where Gaussian masks were applied
    
    % Input validation
    validateattributes(Y, {'numeric'}, {'3d', 'finite', 'nonnan'});
    
    % Get threshold slice
    if nargin < 2 || isempty(threshold_slice)
        threshold_slice = selectSlice(Y);
    end
    validateattributes(threshold_slice, {'numeric'}, {'scalar', 'integer', 'positive', '<=', size(Y,3)});
    
    % Get defect centers
    defect_centers = selectDefectCenters(Y(:,:,threshold_slice));
    if isempty(defect_centers)
        error('No defect centers selected');
    end
    
    % Define radii for each defect
    [sigmas, last_radius] = defineDefectRadii(Y(:,:,threshold_slice), defect_centers);
    
    % Apply Gaussian masks
    [Y_masked, defect_mask] = applyGaussianMasks(Y, defect_centers, sigmas);
    
    % Display results
    displayResults(Y, Y_masked, defect_mask, threshold_slice, defect_centers, sigmas);
end

function slice_num = selectSlice(Y)
    % Let user select a slice for defect selection
    f = figure('Name', 'Select Slice', 'Position', [100, 100, 800, 600]);
    d3gridDisplay(Y, 'dynamic');
    slice_num = input('Enter slice number for selecting defects: ');
    close(f);
end

function centers = selectDefectCenters(slice_data)
    % Let user select defect centers
    f = figure('Name', 'Defect Selection', 'Position', [100, 100, 800, 600]);
    imagesc(slice_data);
    colormap('jet'); colorbar; axis image;
    caxis([min(slice_data(:)), max(slice_data(:))]);
    title('Select defect centers (click points, press Enter when done)');
    
    [x, y] = getpts;
    centers = [x, y];
    close(f);
end

function [sigmas, last_radius] = defineDefectRadii(slice_data, defect_centers)
    % Define radii for each defect center
    sigmas = zeros(size(defect_centers, 1), 1);
    last_radius = [];
    f = figure('Name', 'Define Radii', 'Position', [100, 100, 800, 600]);
    
    for i = 1:size(defect_centers, 1)
        % Display current defect
        clf;
        imagesc(slice_data);
        colormap('jet'); colorbar; axis image;
        caxis([min(slice_data(:)), max(slice_data(:))]);
        
        center = defect_centers(i,:);
        hold on;
        plot(center(1), center(2), 'r+', 'MarkerSize', 15, 'LineWidth', 2);
        plot(center(1), center(2), 'ro', 'MarkerSize', 20, 'LineWidth', 2);
        
        title(sprintf('Define radius for defect %d/%d\nCenter at (%.1f, %.1f)\nDrag circle edge to adjust radius, press Enter when done', ...
            i, size(defect_centers,1), center(1), center(2)));
        
        % Show radius options if not first defect
        if ~isempty(last_radius) && i > 1
            createRadiusButtons(last_radius, center);
        end
        
        % Create and wait for circle
        h = createCircle(center, last_radius);
        if ~isvalid(h)
            error('Circle ROI was deleted before completion');
        end
        
        % Store radius and draw final circle
        radius = h.Radius;
        sigmas(i) = radius/2;
        last_radius = radius;
        
        % Draw final circle
        th = 0:pi/50:2*pi;
        plot(radius * cos(th) + center(1), radius * sin(th) + center(2), 'r--', 'LineWidth', 1.5);
        
        % Handle last defect
        if i == size(defect_centers, 1)
            title('All defects processed!\nPress any key to continue...');
            pause;
        else
            pause(0.5);
        end
    end
    close(f);
end

function createRadiusButtons(last_radius, center)
    % Create buttons for radius options
    use_last_btn = uicontrol('Style', 'pushbutton', ...
        'String', sprintf('Use Last Radius (%.1f)', last_radius), ...
        'Position', [10, 10, 150, 30], ...
        'Callback', @(src,evt) useLastRadius(last_radius, center));
    next_btn = uicontrol('Style', 'pushbutton', ...
        'String', 'Draw New Circle', ...
        'Position', [170, 10, 100, 30], ...
        'Callback', @(src,evt) proceed());
    uiwait;
end

function h = createCircle(center, radius)
    % Create a circle ROI
    if ~isempty(radius)
        h = drawcircle('Center', center, 'Radius', radius, 'Color', 'r', 'FaceAlpha', 0.1);
    else
        h = drawcircle('Center', center, 'Color', 'r', 'FaceAlpha', 0.1);
    end
    addlistener(h, 'ROIMoved', @(src,evt) forceCenter(src, center));
    wait(h);
end

function [Y_masked, defect_mask] = applyGaussianMasks(Y, defect_centers, sigmas)
    % Apply Gaussian masks to the data
    Y_masked = Y;
    defect_mask = false(size(Y));
    [X, Y_coord] = meshgrid(1:size(Y,2), 1:size(Y,1));
    
    for i = 1:size(Y, 3)
        mask = ones(size(Y,1), size(Y,2));
        for j = 1:size(defect_centers, 1)
            r0 = defect_centers(j,:);
            gaussian = 0.99 * exp(-((X-r0(1)).^2 + (Y_coord-r0(2)).^2)/(2*sigmas(j)^2));
            mask = mask .* (1 - gaussian);
        end
        Y_masked(:,:,i) = Y(:,:,i) .* mask;
        defect_mask(:,:,i) = mask < 1;
    end
end

function displayResults(Y, Y_masked, defect_mask, threshold_slice, defect_centers, sigmas)
    % Display the results
    figure('Name', 'Masking Result', 'Position', [100, 100, 1200, 500]);
    
    subplot(1,3,1);
    imagesc(Y(:,:,threshold_slice));
    colorbar; title('Original Data'); axis image;
    caxis([min(Y(:)), max(Y(:))]);
    
    subplot(1,3,2);
    imagesc(1 - defect_mask(:,:,threshold_slice));
    colorbar; title('Gaussian Mask'); axis image;
    colormap(gca, 'gray');
    
    subplot(1,3,3);
    imagesc(Y_masked(:,:,threshold_slice));
    colorbar; title('After Gaussian Masking'); axis image;
    caxis([min(Y_masked(:)), max(Y_masked(:))]);
    colormap(gca, 'jet');
    
    % Print parameters
    fprintf('\nGaussian Masking Parameters:\n');
    fprintf('Number of defects masked: %d\n', size(defect_centers, 1));
    for i = 1:size(defect_centers, 1)
        fprintf('Defect %d: position (%.1f, %.1f), σ = %.1f pixels\n', ...
            i, defect_centers(i,1), defect_centers(i,2), sigmas(i));
    end
end

function useLastRadius(radius, center)
    % Handle using last radius
    existing_circles = findobj(gca, 'Type', 'images.roi.Circle');
    if ~isempty(existing_circles)
        delete(existing_circles);
    end
    h = createCircle(center, radius);
    if isvalid(h)
        uiresume;
    else
        error('Circle ROI was deleted before completion');
    end
end

function proceed()
    % Handle proceeding to next defect
    uiresume;
end

function forceCenter(h, center)
    % Force circle to stay centered
    if isvalid(h)
        h.Center = center;
    end
end 