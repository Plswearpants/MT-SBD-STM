function [quality_metric, residual] = computeResidualQuality(Y_total, reconstructed_kernels, reconstructed_activations, noise_var)
    % Calculate total reconstruction
    Y_reconstructed = zeros(size(Y_total));
    for k = 1:length(reconstructed_kernels)
        Y_reconstructed = Y_reconstructed + convfft2(reconstructed_kernels{k}, reconstructed_activations(:,:,k));
    end
    
    % Normalize both signals
    Y_total_norm = proj2oblique(Y_total);
    Y_reconstructed_norm = proj2oblique(Y_reconstructed);
    
    % Calculate residual using normalized signals
    residual = Y_total_norm - Y_reconstructed_norm;
    residual_var = var(residual, 0, 'all');
    
    % Calculate quality metric
    quality_metric = noise_var / residual_var;
    
    % Optional: Print diagnostic information
    fprintf('var(noise) = %.4e\n', noise_var);
    fprintf('var(residual) = %.4e\n', residual_var);
    fprintf('var(noise)/var(residual) = %.4f\n', quality_metric);
end 