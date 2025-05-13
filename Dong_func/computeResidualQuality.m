function [quality_metric, residual] = computeResidualQuality(Y_total, reconstructed_kernels, reconstructed_activations, noise_var)
    % Calculate total reconstruction
    Y_reconstructed = zeros(size(Y_total));
    if size(Y_total,3) == 1
        for k = 1:length(reconstructed_kernels)
            Y_reconstructed = Y_reconstructed + convfft2(reconstructed_kernels{k}, reconstructed_activations(:,:,k));
        end

        
        % Calculate residual using normalized signals
        residual = Y_total - Y_reconstructed;
        residual_var = var(residual, 0, 'all');
        
        % Calculate quality metric
        quality_metric = noise_var / residual_var;
        
        % Optional: Print diagnostic information
        fprintf('var(noise) = %.4e\n', noise_var);
        fprintf('var(residual) = %.4e\n', residual_var);
        fprintf('var(noise)/var(residual) = %.4f\n', quality_metric);
    else
        for k = 1:length(reconstructed_kernels)
            X_rec = reconstructed_activations(:,:,k);
            Y_reconstructed = Y_reconstructed + convfft3(reconstructed_kernels{k}, X_rec(:,:,ones(1,size(Y_total,3))));
        end

        % Calculate residual using normalized signals
        residual = Y_total - Y_reconstructed;
        residual_var = reshape(var(residual, 0, [1,2]),[1,size(Y_total,3)]);
        
        % Calculate quality metric
        quality_metric = mean(noise_var ./ residual_var);
    end

end 