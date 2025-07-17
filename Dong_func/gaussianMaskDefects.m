function [Y_masked, combined_mask, defect_centers, sigmas] = gaussianMaskDefects(Y, threshold_slice, num_defect_type, varargin)
%GAUSSIANMASKDEFECTS Mask defects in 3D data using Gaussian masks.
%   [Y_masked, combined_mask, defect_centers, sigmas] = gaussianMaskDefects(Y, threshold_slice, num_defect_type)
%   [Y_masked, combined_mask, defect_centers, sigmas] = gaussianMaskDefects(Y, threshold_slice, num_defect_type, mask)
%   - If mask is provided, it is used directly. Otherwise, mask is created interactively.
%   - combined_mask is a 2D mask (height x width) applied to all slices.

    validateattributes(Y, {'numeric'}, {'3d', 'finite', 'nonnan'});
    if nargin < 2 || isempty(threshold_slice)
        threshold_slice = selectSlice(Y);
    end
    validateattributes(threshold_slice, {'numeric'}, {'scalar', 'integer', 'positive', '<=', size(Y,3)});
    if nargin >= 4 && ~isempty(varargin{1})
        combined_mask = varargin{1};
        defect_centers = {};
        sigmas = [];
    else
        Y_temp = Y;
        defect_centers = cell([1,num_defect_type]);
        sigmas = cell([1,num_defect_type]);
        combined_mask = ones(size(Y,1), size(Y,2));
        for i = 1:num_defect_type
            fprintf('Select activation for defect type %d\n',i);
            defect_centers{i} = selectDefectCenters(Y_temp(:,:,threshold_slice));
            if isempty(defect_centers{i})
                disp('No defect centers selected');
                continue;
            end
            sigmas{i} = defineDefectRadii(Y(:,:,threshold_slice), defect_centers{i});
            % Build mask for this defect type
            mask_i = buildGaussianMask(size(Y,1), size(Y,2), defect_centers{i}, sigmas{i});
            combined_mask = combined_mask .* mask_i;
            % Display result after each type
            for k = 1:size(Y,3)
                Y_temp(:,:,k) = Y(:,:,k) .* combined_mask;
            end
            displayResults(Y, Y_temp, threshold_slice, defect_centers{i}, sigmas{i});
        end
    end
    % Apply the mask to all slices
    Y_masked = Y;
    for i = 1:size(Y,3)
        Y_masked(:,:,i) = Y(:,:,i) .* combined_mask;
    end
end

function mask = buildGaussianMask(h, w, centers, sigmas)
    mask = ones(h, w);
    [X, Y_coord] = meshgrid(1:w, 1:h);
    for j = 1:size(centers, 1)
        r0 = centers(j,:);
        distance_squared = (X-r0(1)).^2 + (Y_coord-r0(2)).^2;
        distance = sqrt(distance_squared);
        sigma = sigmas(j);
        step_loc = 2*sigma;
        step_shapeness = 10;
        smooth_step_fn = @(d, loc, shape) 0.5 + 0.5*(tanh(-shape*(d-loc)));
        smooth_step = smooth_step_fn(distance, step_loc, step_shapeness);
        gaussian = 0.99 * exp(-distance_squared/(2*sigma^2)) .* smooth_step;
        mask = mask .* (1 - gaussian);
    end
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
    colormap('gray'); colorbar; axis image;
    clim([min(slice_data(:)), max(slice_data(:))]);
    title('Select defect centers (click points, press Enter when done)');
    
    [x, y] = getpts;
    centers = [x, y];
    close(f);
end

function sigmas = defineDefectRadii(slice_data, defect_centers)
    % Define radii for each defect center
    sigmas = zeros(size(defect_centers, 1), 1);
    f = figure('Name', 'Define Radii', 'Position', [100, 100, 800, 600]);
    
    % Default radius (can be adjusted as needed)
    default_radius = 0.1;
    pre_rad = default_radius;
    for i = 1:size(defect_centers, 1)
        % Display current defect
        clf;
        imagesc(slice_data);
        colormap('gray'); colorbar; axis image;
        caxis([min(slice_data(:)), max(slice_data(:))]);
        
        center = defect_centers(i,:);
        hold on;
        plot(center(1), center(2), 'r+', 'MarkerSize', 1, 'LineWidth', 0.2);
        
        % Create radius display and edit box
        radius_text = uicontrol('Style', 'text', ...
            'String', 'Radius: ', ...
            'Position', [10, 10, 60, 20]);
        radius_edit = uicontrol('Style', 'edit', ...
            'String', num2str(pre_rad), ...
            'Position', [70, 10, 60, 20], ...
            'Callback', @(src,evt) updateRadius(src, center));
        
        % Create Confirm button
        confirm_button = uicontrol('Style', 'pushbutton', ...
            'String', 'Confirm', ...
            'Position', [140, 10, 80, 20], ...
            'Callback', @(src,evt) uiresume(f));
        
        % Create circle with default radius
        h = drawcircle('Center', center, 'Radius', pre_rad, 'Color', 'r', 'FaceAlpha', 0.1);
        
        % Add listeners to keep circle centered and update radius display
        addlistener(h, 'ROIMoved', @(src,evt) handleROIMoved(src, center, radius_edit));
        
        % Wait for button press
        uiwait(f);
        
        % Get final radius and compute sigma
        radius = h.Radius;
        sigmas(i) = radius;
        pre_rad = radius;
        % Draw final circle
        th = 0:pi/50:2*pi;
        plot(radius * cos(th) + center(1), radius * sin(th) + center(2), 'r--', 'LineWidth', 1.5);
    end
    close(f);
    
    function handleROIMoved(src, center, radius_edit)
        % Handle both center and radius changes
        if isvalid(src)
            % Force center to stay at defect point
            src.Center = center;
            % Update radius display
            radius_edit.String = num2str(src.Radius, '%.1f');
        end
    end
    
    function updateRadius(src, center)
        % Update circle when radius is edited
        new_radius = str2double(src.String);
        if ~isnan(new_radius) && new_radius > 0
            delete(findobj(gca, 'Type', 'images.roi.Circle'));
            h = drawcircle('Center', center, 'Radius', new_radius, 'Color', 'r', 'FaceAlpha', 0.1);
            addlistener(h, 'ROIMoved', @(src,evt) handleROIMoved(src, center, radius_edit));
        end
    end
end

function displayResults(Y, Y_masked, threshold_slice, defect_centers, sigmas)
    invgray = flipud(gray);
    % Display the results
    figure('Name', 'Masking Result', 'Position', [100, 100, 1800, 800]);
    Y_slice = Y(:,:,threshold_slice);
    Y_masked_slice = Y_masked(:,:,threshold_slice);
    Y_inverse_mask_slice = Y_slice - Y_masked_slice;

    % Calculate QPI (2D FFT)
    Y_fft = fftshift(abs(fft2(Y_slice)));
    Y_masked_fft = fftshift(abs(fft2(Y_masked_slice)));
    Y_inverse_fft = fftshift(abs(fft2(Y_inverse_mask_slice)));
    
    % Normalize FFT for better visualization
    Y_fft = log10(Y_fft + 1);
    Y_masked_fft = log10(Y_masked_fft + 1);
    Y_inverse_fft = log10(Y_inverse_fft + 1);

    % Real space images (first row)
    subplot(2,3,1);
    imagesc(Y_slice);
    colorbar; title('Original Data'); axis image;
    caxis([min(Y_slice(:)), max(Y_slice(:))]);
    colormap(gca, invgray);
    
    subplot(2,3,2);
    imagesc(Y_masked_slice);
    colorbar; title('After Gaussian Masking'); axis image;
    caxis([min(Y_slice(:)), max(Y_slice(:))]);
    colormap(gca, invgray);
    
    subplot(2,3,3);
    imagesc(Y_inverse_mask_slice);
    colorbar; title('Inverse mask'); axis image;
    caxis([min(Y_slice(:)), max(Y_slice(:))]);
    colormap(gca, invgray);
    
    % QPI images (second row)
    subplot(2,3,4);
    imagesc(Y_fft);
    colorbar; title('QPI (Original)'); axis image;
    colormap(gca, invgray);
    caxis auto;
    
    subplot(2,3,5);
    imagesc(Y_masked_fft);
    colorbar; title('QPI (Masked)'); axis image;
    colormap(gca, invgray);
    caxis auto;
    
    subplot(2,3,6);
    imagesc(Y_inverse_fft);
    colorbar; title('QPI (Inverse)'); axis image;
    colormap(gca, invgray);
    caxis auto;

    % Print parameters
    fprintf('\nGaussian Masking Parameters:\n');
    fprintf('Number of defects masked: %d\n', size(defect_centers, 1));
    for i = 1:size(defect_centers, 1)
        fprintf('Defect %d: position (%.1f, %.1f), σ = %.1f pixels\n', ...
            i, defect_centers(i,1), defect_centers(i,2), sigmas(i));
    end
end 