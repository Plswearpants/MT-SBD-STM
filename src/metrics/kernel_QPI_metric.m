function [similarity_score, qpi_output, qpi_gt] = kernel_QPI_metric(output_kernel, gt_kernel)
    % Computes similarity between output kernel and ground truth kernel in QPI space
    %
    % Inputs:
    %   output_kernel: The output kernel from the algorithm
    %   gt_kernel: Ground truth kernel
    %
    % Outputs:
    %   similarity_score: Metric value indicating similarity (higher is better)
    %   qpi_output: QPI of output kernel (for visualization)
    %   qpi_gt: QPI of ground truth kernel (for visualization)
    
    % Convert both kernels to QPI space
    qpi_output = abs(fftshift(fft2(output_kernel)));
    qpi_gt = abs(fftshift(fft2(gt_kernel)));
    
    % Normalize QPIs to account for potential intensity differences
    qpi_output = qpi_output / max(qpi_output(:));
    qpi_gt = qpi_gt / max(qpi_gt(:));
    
    % Calculate similarity score
    % Using normalized cross-correlation (values between -1 and 1, 1 being perfect match)
    similarity_score = corr2(qpi_output, qpi_gt);
    
    % Alternative metrics could be:
    % MSE = mean((qpi_output(:) - qpi_gt(:)).^2);
    % SSIM = ssim(qpi_output, qpi_gt);
end 