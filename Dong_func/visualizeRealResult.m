function visualizeRealResult(Y_ref, A_ref, X_ref, b_ref, extras_ref)
    % Visualize results from SBD-STM reconstruction for real data
    %
    % Inputs:
    %   A_ref: Reconstructed kernels (cell array)
    %   X_ref: Reconstructed activation maps
    %   b_ref: Bias terms [num_kernels x 1]
    %   extras_ref: Struct containing residuals and other metrics
    
    num_kernels = length(A_ref);
    
    % Get original observation and residual
    Y = Y_ref;  % Original observation
    Y_res = extras_ref.phase1.residuals(:,:,end);  % Final residual
    
    % Compute reconstructed image
    Y_rec = zeros(size(Y));
    for i = 1:num_kernels
        Y_rec = Y_rec + convfft2(A_ref{i}, X_ref(:,:,i));
        Y_rec = Y_rec + b_ref(i);  % Add bias term if needed
    end
    
    % 1. Original vs Reconstructed vs Residual (1x3 plot)
    figure('Name', 'Image Comparison');
    subplot(131);
    imagesc(Y);
    title('Original Image');
    colorbar;
    axis image;
    
    subplot(132);
    imagesc(Y_rec);
    title('Reconstructed Image');
    colorbar;
    axis image;
    
    subplot(133);
    imagesc(Y_res);
    title('Residual Image');
    colorbar;
    axis image;
    
    % Set same color limits for all three plots
    clim = [min(Y(:)), max(Y(:))];
    for i = 1:3
        subplot(1,3,i);
        caxis(clim);
    end
    sgtitle('Original vs Reconstructed vs Residual');

    % print the quality metric = var(noise)/var(residual)  
    quality_metric = extras_ref.phase1.quality_metrics(end);
    fprintf('Quality metric = var(noise)/var(residual) = %f\n', quality_metric);
    
    % 2. Auto-contrasted Residual
    figure('Name', 'Auto-contrasted Residual');
    imagesc(Y_res);
    title('Auto-contrasted Residual');
    colorbar;
    axis image;
    
    % 3. Q-space Analysis
    figure('Name', 'Q-space Analysis');
    
    % Calculate Q-space for original image
    Q_Y = qpiCalculate(Y);
    Q_Y_avg = mean(Q_Y, 3);  % Average over voltage points if 3D
    
    % Plot Q-space of original image
    subplot(2, num_kernels + 1, num_kernels + 2);
    imagesc(Q_Y_avg);
    title('Q-space (Original)');
    colorbar;
    axis image;
    
    % Calculate and plot Q-space for each kernel
    for k = 1:num_kernels
        % Calculate Q-space for kernel
        Q_kernel = qpiCalculate(A_ref{k});
        Q_kernel_avg = mean(Q_kernel, 3);  % Average over voltage points if 3D
    
        subplot(2, num_kernels + 1, k + 1);
        imagesc(A_ref{k});
        title(['Kernel ' num2str(k) ' (Original)']);
        colorbar;
        axis square;

        subplot(2, num_kernels + 1, num_kernels+ k + 2);
        imagesc(Q_kernel_avg);
        title(['Q-space (Kernel ' num2str(k) ')']);
        colorbar;
        axis image;
    end
    sgtitle('Q-space Analysis');
end 