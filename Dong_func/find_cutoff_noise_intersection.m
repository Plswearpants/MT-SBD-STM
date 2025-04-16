function [M, M_pixels,fig_handle] = find_cutoff_noise_intersection(rho_single, SNR, N, show_plots)
    % Determines cutoff M by finding where signal meets noise level in radial directions
    % 
    % Inputs:
    %   rho_single: 2D array of QPI pattern
    %   SNR: Signal-to-noise ratio
    %   N: Lattice span in one direction (NxN lattice)
    %   show_plots: boolean to control visualization (default true)
    %
    % Outputs:
    %   M: Cutoff radius in lattice units
    %   M_pixels: Cutoff radius in pixels
    %   fig_handle: Handle to the figure (if plots shown)
    
    if nargin < 4
        show_plots = true;
    end
    
    % Get pixel to lattice conversion factor
    [ny, nx] = size(rho_single);
    pixels_per_lattice = nx/N; % Assuming square image and lattice
    
    % Calculate noise level based on signal variance
    signal_variance = var(rho_single(:));
    noise_level = sqrt(signal_variance/SNR);
    
    % Get center of the pattern
    center = ceil([ny, nx]/2);
    
    % Calculate radial profile
    [r, avg_intensity] = get_radial_profile(rho_single, center);
    signal_envelope = abs(avg_intensity);
    
    % Convert r to lattice units for peak finding
    r_lattice = r/pixels_per_lattice;
    
    % Find peaks in the signal
    [peak_values, peak_locs] = findpeaks(signal_envelope, ...
        'MinPeakProminence', noise_level/5); % Minimum peak prominence
    
    % Find first peak that falls below half of noise level
    peak_below_noise = find(peak_values < noise_level/2, 1, 'first');
    
    if isempty(peak_below_noise)
        M_pixels = min(nx, ny)/2;
    else
        M_pixels = peak_locs(peak_below_noise);
    end
    
    % Convert M to lattice units
    M = M_pixels/pixels_per_lattice;
    
    % Ensure M is reasonable in lattice units
    min_size_lattice = 2;
    max_size_lattice = floor(N/2);
    M = max(min_size_lattice, min(M, max_size_lattice));
    
    % Add noise to data for visualization
    noise = sqrt(signal_variance/SNR) * randn(size(rho_single));
    noisy_data = rho_single + noise;
    
    if show_plots
        fig_handle = figure('Position', [100 100 1500 400]);
        
        % Plot 1: Original QPI pattern
        subplot(1,4,1)
        imagesc(rho_single);
        colorbar;
        title('Original QPI Pattern');
        axis square
        hold on
        % Draw circle at cutoff radius (in pixels)
        th = 0:pi/50:2*pi;
        xunit = M*pixels_per_lattice * cos(th) + center(2);
        yunit = M*pixels_per_lattice * sin(th) + center(1);
        plot(xunit, yunit, 'r--', 'LineWidth', 2);
        hold off
        
        % Plot 2: Noisy QPI pattern
        subplot(1,4,2)
        imagesc(noisy_data);
        colorbar;
        title(sprintf('Noisy QPI Pattern (SNR=%.1f)', SNR));
        axis square
        hold on
        % Draw same cutoff circle
        plot(xunit, yunit, 'r--', 'LineWidth', 2);
        hold off
        
        % Plot 3: Radial profile with noise level and peaks
        subplot(1,4,3)
        plot(r_lattice, avg_intensity, 'b-', 'LineWidth', 1, 'DisplayName', 'Signal');
        hold on
        plot(r_lattice, signal_envelope, 'r-', 'LineWidth', 1, 'DisplayName', 'Signal Envelope');
        % Plot all detected peaks
        plot(peak_locs/pixels_per_lattice, peak_values, 'ko', 'MarkerFaceColor', 'y', ...
            'DisplayName', 'Peaks');
        % Plot noise levels
        plot(xlim, [noise_level noise_level], 'g--', 'LineWidth', 1.5, ...
            'DisplayName', 'Noise Level');
        plot(xlim, [-noise_level -noise_level], 'g--', 'LineWidth', 1.5, ...
            'HandleVisibility', 'off');
        % Plot cutoff
        plot([M M], ylim, 'k--', 'LineWidth', 1.5, 'DisplayName', 'Cutoff');
        title('Radial Profile with Peaks and Noise Level');
        xlabel('Radius (lattice units)');
        ylabel('Intensity');
        grid on;
        legend('Location', 'best');
        
        % Plot 4: Cropped QPI pattern
        subplot(1,4,4)
        range_y = max(1, center(1)-M*pixels_per_lattice):min(ny, center(1)+M*pixels_per_lattice);
        range_x = max(1, center(2)-M*pixels_per_lattice):min(nx, center(2)+M*pixels_per_lattice);
        range_y = round(range_y);
        range_x = round(range_x);
        cropped = rho_single(range_y, range_x);
        imagesc(cropped);
        colorbar;
        title(sprintf('Cropped Pattern (M=%.1f lattice units)', M));
        axis square
        
        sgtitle(sprintf('QPI Pattern Analysis (SNR=%.1f)', SNR));
    else
        fig_handle = [];
    end
end

function [r, avg_intensity] = get_radial_profile(img, center)
    [X, Y] = meshgrid(1:size(img,2), 1:size(img,1));
    R = sqrt((X-center(2)).^2 + (Y-center(1)).^2);
    R = round(R);
    
    max_r = max(R(:));
    r = 0:max_r;
    avg_intensity = zeros(size(r));
    
    for i = 1:length(r)
        mask = (R == r(i));
        if any(mask(:))
            avg_intensity(i) = mean(img(mask));
        end
    end
    
    % Smooth the profile slightly to reduce noise but preserve peaks
    avg_intensity = smoothdata(avg_intensity, 'gaussian', 3);
end
