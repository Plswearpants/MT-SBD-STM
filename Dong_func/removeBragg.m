function [Y_removed, mask2d, recipe] = removeBragg(Y, opts)
%REMOVEBRAGG Remove Bragg peaks in reciprocal space (interactive or replay).
%   Y_removed = removeBragg(Y)
%   [Y_removed, mask2d, recipe] = removeBragg(Y)
%   Y_removed = removeBragg(Y, opts)
%   [Y_removed, mask2d, recipe] = removeBragg(Y, opts)
%
%   When opts is a struct with non-empty opts.mask2d, runs in non-interactive
%   replay mode: applies opts.mask2d to FFT of Y and returns Y_removed. No UI.
%   When opts is absent or opts.mask2d is empty, runs interactively and
%   optionally returns mask2d and recipe for checkpoint recording.

    if nargin < 2
        opts = struct();
    end
    if ~isstruct(opts)
        opts = struct();
    end

    % Non-interactive replay: apply provided mask
    if isfield(opts, 'mask2d') && ~isempty(opts.mask2d)
        mask2d = opts.mask2d;
        recipe = [];
        if isfield(opts, 'recipe')
            recipe = opts.recipe;
        end
        [rows, cols, nsl] = size(Y);
        if ~isequal(size(mask2d), [rows, cols])
            error('removeBragg: opts.mask2d size [%d,%d] does not match data [%d,%d].', ...
                size(mask2d,1), size(mask2d,2), rows, cols);
        end
        QPI = zeros(size(Y));
        for i = 1:nsl
            QPI(:,:,i) = fftshift(fft2(Y(:,:,i)));
        end
        Y_removed = Y;
        for i = 1:nsl
            QPI(:,:,i) = QPI(:,:,i) .* mask2d;
            Y_removed(:,:,i) = real(ifft2(ifftshift(QPI(:,:,i))));
        end
        if nargout < 2
            clear mask2d recipe;
        end
        return;
    end

    % Interactive path
    QPI = zeros(size(Y));
    for i = 1:size(Y,3)
        QPI(:,:,i) = fftshift(fft2(Y(:,:,i)));
    end
    QPI_logabs = (abs(QPI));

    figure;
    d3gridDisplay(QPI_logabs, 'dynamic');
    slice = input('Enter the slice number to remove the bragg peaks: ');
    close;

    figure;
    imagesc(QPI_logabs(:,:,slice)); axis image; colormap hot; title('Select Bragg peaks to apply Gaussian window');

    num_peaks = input('Enter the number of unique Bragg peaks to process: ');
    Y_removed = Y;
    QPI_removed = QPI;

    [rows, cols, ~] = size(QPI);
    center_row = floor(rows/2)+1;
    center_col = floor(cols/2)+1;

    mask = ones(rows, cols);

    removal_method = input('Choose removal method (1 for Gaussian window, 2 for complete removal): ');

    for i = 1:num_peaks
        disp(['Select Bragg peak #', num2str(i), ' with elliptical ROI']);

        h = drawellipse('Color', 'b');
        wait(h);

        center = h.Center;
        semi_axes = h.SemiAxes;
        rotation = h.RotationAngle;

        if length(semi_axes) ~= 2
            error('Semi-axes must be a 2-element vector for independent control');
        end

        disp(['Current semi-axes: [', num2str(semi_axes(1)), ', ', num2str(semi_axes(2)), ']']);
        disp('You can adjust these values independently using the ellipse handles');

        sym_center_x = 2*center_row - center(1);
        sym_center_y = 2*center_col - center(2);

        hold on;
        h_marker = plot(center(1), center(2), 'g+', 'MarkerSize', 10);
        sym_marker = plot(sym_center_x, sym_center_y, 'g+', 'MarkerSize', 10);

        sym_h = drawellipse('Center', [sym_center_x, sym_center_y], ...
                           'SemiAxes', semi_axes, ...
                           'RotationAngle', rotation, ...
                           'Color', 'b');

        h.InteractionsAllowed = 'all';
        sym_h.InteractionsAllowed = 'all';

        addlistener(h, 'MovingROI', @(src,evt) updateEllipseAndMarker(src, sym_h, sym_marker, h_marker, rows, cols));
        addlistener(sym_h, 'MovingROI', @(src,evt) updateEllipseAndMarker(src, h, h_marker, sym_marker, rows, cols));

        disp('Adjust the ellipses as needed. Press any key when done...');
        pause;

        center = h.Center;
        semi_axes = h.SemiAxes;
        rotation = h.RotationAngle;
        sym_center = sym_h.Center;

        [X_grid, Y_grid] = meshgrid(1:cols, 1:rows);

        theta = deg2rad(rotation);
        X_rot = (X_grid - center(1)) * cos(theta) + (Y_grid - center(2)) * sin(theta);
        Y_rot = -(X_grid - center(1)) * sin(theta) + (Y_grid - center(2)) * cos(theta);

        ellipse_mask = (X_rot.^2/semi_axes(1)^2 + Y_rot.^2/semi_axes(2)^2) <= 1;

        if removal_method == 1
            sigma_x = semi_axes(1);
            sigma_y = semi_axes(2);
            gaussian_window = exp(-(X_rot.^2/(2*sigma_x^2) + Y_rot.^2/(2*sigma_y^2)));
            peak_mask = 1 - gaussian_window .* ellipse_mask;
        else
            peak_mask = 1 - ellipse_mask;
        end

        X_sym_rot = (X_grid - sym_center(1)) * cos(theta) + (Y_grid - sym_center(2)) * sin(theta);
        Y_sym_rot = -(X_grid - sym_center(1)) * sin(theta) + (Y_grid - sym_center(2)) * cos(theta);

        if removal_method == 1
            sym_gaussian_window = exp(-(X_sym_rot.^2/(2*sigma_x^2) + Y_sym_rot.^2/(2*sigma_y^2)));
            sym_ellipse_mask = (X_sym_rot.^2/semi_axes(1)^2 + Y_sym_rot.^2/semi_axes(2)^2) <= 1;
            sym_peak_mask = 1 - sym_gaussian_window .* sym_ellipse_mask;
        else
            sym_ellipse_mask = (X_sym_rot.^2/semi_axes(1)^2 + Y_sym_rot.^2/semi_axes(2)^2) <= 1;
            sym_peak_mask = 1 - sym_ellipse_mask;
        end

        mask = mask .* peak_mask .* sym_peak_mask;

        delete(h);
        delete(sym_h);
        delete(h_marker);
        delete(sym_marker);
    end
    close;

    for i = 1:size(QPI, 3)
        QPI_removed(:,:,i) = QPI(:,:,i) .* mask;
        Y_removed(:,:,i) = real(ifft2(ifftshift(QPI_removed(:,:,i))));
    end

    % Provenance for checkpoint
    recipe = struct('slice', slice, 'num_peaks', num_peaks, 'removal_method', removal_method);
    mask2d = mask;

    figure;
    subplot(2,2,1); imagesc(Y(:,:,slice)); axis image; title('Original Image');
    subplot(2,2,2); imagesc(log(abs(QPI(:,:,slice)))); axis image; title('Original FFT');
    subplot(2,2,3); imagesc(Y_removed(:,:,slice)); axis image; title('Filtered Image');
    subplot(2,2,4); imagesc(log(abs(QPI_removed(:,:,slice)))); axis image; title('Filtered FFT with Gaussian Window');

    if nargout < 2
        clear mask2d recipe;
    end
end

function updateEllipseAndMarker(src, target, marker, src_marker, rows, cols)
    dist_to_left = src.Center(1) - 1;
    dist_to_top = src.Center(2) - 1;
    sym_center = [cols - dist_to_left, rows - dist_to_top];
    target.Center = sym_center;
    target.SemiAxes = src.SemiAxes;
    target.RotationAngle = src.RotationAngle;
    src_marker.XData = src.Center(1);
    src_marker.YData = src.Center(2);
    marker.XData = sym_center(1);
    marker.YData = sym_center(2);
end
