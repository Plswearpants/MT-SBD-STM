function metrics = compareSharpness(original, reconstructed)
    % Compare sharpness between original and reconstructed signals
    % 
    % Inputs:
    %   original: Original signal/kernel
    %   reconstructed: Reconstructed signal/kernel
    %
    % Outputs:
    %   metrics: Struct containing different sharpness metrics
    %     .gradient_magnitude: Average gradient magnitude
    %     .laplacian_variance: Variance of Laplacian
    %     .frequency_content: High-frequency content ratio
    
    metrics = struct();
    
    % 1. Gradient-based sharpness
    [Gx_orig, Gy_orig] = gradient(original);
    [Gx_recon, Gy_recon] = gradient(reconstructed);
    
    metrics.gradient_magnitude.original = mean(sqrt(Gx_orig.^2 + Gy_orig.^2), 'all');
    metrics.gradient_magnitude.reconstructed = mean(sqrt(Gx_recon.^2 + Gy_recon.^2), 'all');
    metrics.gradient_magnitude.ratio = metrics.gradient_magnitude.reconstructed / ...
                                     metrics.gradient_magnitude.original;
    
    % 2. Laplacian-based sharpness
    lap_orig = del2(original);
    lap_recon = del2(reconstructed);
    
    metrics.laplacian_variance.original = var(lap_orig, 0, 'all');
    metrics.laplacian_variance.reconstructed = var(lap_recon, 0, 'all');
    metrics.laplacian_variance.ratio = metrics.laplacian_variance.reconstructed / ...
                                     metrics.laplacian_variance.original;
    
    % 3. Frequency domain analysis
    F_orig = fftshift(fft2(original));
    F_recon = fftshift(fft2(reconstructed));
    
    % Create frequency mask (consider frequencies above 50% as high)
    [M, N] = size(original);
    [X, Y] = meshgrid(-N/2:N/2-1, -M/2:M/2-1);
    R = sqrt(X.^2 + Y.^2);
    high_freq_mask = R > (min(M,N)/4);  % Adjust threshold as needed
    
    metrics.frequency_content.original = sum(abs(F_orig(high_freq_mask)), 'all') / ...
                                       sum(abs(F_orig(:)), 'all');
    metrics.frequency_content.reconstructed = sum(abs(F_recon(high_freq_mask)), 'all') / ...
                                            sum(abs(F_recon(:)), 'all');
    metrics.frequency_content.ratio = metrics.frequency_content.reconstructed / ...
                                    metrics.frequency_content.original;
end 