function metrics = compareNoise(kernel1, kernel2)
    % Compare noise levels between two kernels using frequency analysis
    % 
    % Inputs:
    %   kernel1, kernel2: Two kernels to compare (real space)
    %
    % Outputs:
    %   metrics: Struct containing different noise metrics
    
    % Input validation
    validateattributes(kernel1, {'numeric'}, {'2d', 'real', 'finite'});
    validateattributes(kernel2, {'numeric'}, {'2d', 'real', 'finite'});
    if ~isequal(size(kernel1), size(kernel2))
        error('Kernels must have the same size');
    end
    
    % 1. Basic variance in real space
    kernel1_norm = kernel1 / norm(kernel1(:));
    kernel2_norm = kernel2 / norm(kernel2(:));
    
    metrics.variance.kernel1 = var(kernel1_norm(:));
    metrics.variance.kernel2 = var(kernel2_norm(:));
    metrics.variance.ratio = metrics.variance.kernel2 / metrics.variance.kernel1;
    
    % 2. Frequency domain analysis
    F1 = fftshift(fft2(kernel1_norm));
    F2 = fftshift(fft2(kernel2_norm));
    
    % Plot FFT and let user select signal region
    [M, N] = size(kernel1);
    center_x = N/2;
    center_y = M/2;
    
    figure('Name', 'Adjust Cutoff Frequency');
    imagesc(log(abs(F1) + 1));
    colormap('jet');
    colorbar;
    title('Adjust Square Size to Set Cutoff Frequency');
    axis image;
    
    % Create centered square ROI
    initial_size = min(M,N)/4;  % Start with 1/4 of image size
    pos = [center_x - initial_size/2, center_y - initial_size/2, initial_size, initial_size];
    
    h_rect = drawrectangle('Position', pos, ...
                          'Label', 'Signal Region', ...
                          'Color', 'r', ...
                          'FixedAspectRatio', true, ...  % Force square shape
                          'Rotatable', false);           % Prevent rotation
    
    % Add constraint to keep square centered
    addlistener(h_rect, 'ROIMoved', @(src,evt) centerROI(src, center_x, center_y));
    
    wait(h_rect);
    
    % Get final position and create masks
    final_pos = h_rect.Position;
    cutoff_freq = final_pos(3)/2;  % Half of square size
    
    % Create radial distance matrix
    [X, Y] = meshgrid(-N/2:N/2-1, -M/2:M/2-1);
    R = sqrt(X.^2 + Y.^2);
    
    signal_mask = R <= cutoff_freq;
    noise_mask = ~signal_mask;
    close;
    
    % Power spectrum calculation
    P1 = abs(F1).^2;
    P2 = abs(F2).^2;
    
    % Noise power in high frequencies
    metrics.noise_power.kernel1 = sum(P1(noise_mask), 'all') / sum(P1, 'all');
    metrics.noise_power.kernel2 = sum(P2(noise_mask), 'all') / sum(P2, 'all');
    metrics.noise_power.ratio = metrics.noise_power.kernel2 / metrics.noise_power.kernel1;
    
    % 3. Local variance analysis
    window_size = min(3, min(size(kernel1))-1);
    kernel1_local_var = nlfilter(kernel1_norm, [window_size window_size], @(x) var(x(:)));
    kernel2_local_var = nlfilter(kernel2_norm, [window_size window_size], @(x) var(x(:)));
    
    metrics.local_variance.kernel1 = mean(kernel1_local_var, 'all');
    metrics.local_variance.kernel2 = mean(kernel2_local_var, 'all');
    metrics.local_variance.ratio = metrics.local_variance.kernel2 / metrics.local_variance.kernel1;
    
    % 4. SNR estimation
    metrics.snr.kernel1 = sum(P1(signal_mask), 'all') / sum(P1(noise_mask), 'all');
    metrics.snr.kernel2 = sum(P2(signal_mask), 'all') / sum(P2(noise_mask), 'all');
    metrics.snr.ratio = metrics.snr.kernel2 / metrics.snr.kernel1;
    
    % Store cutoff information
    metrics.cutoff.frequency = cutoff_freq;
    metrics.cutoff.relative = cutoff_freq/(min(M,N)/2);  % Relative to Nyquist
    metrics.cutoff.signal_mask = signal_mask;
    metrics.cutoff.noise_mask = noise_mask;
end

function centerROI(roi, center_x, center_y)
    % Helper function to keep ROI centered
    pos = roi.Position;
    size = pos(3);  % Square size
    new_pos = [center_x - size/2, center_y - size/2, size, size];
    roi.Position = new_pos;
end 