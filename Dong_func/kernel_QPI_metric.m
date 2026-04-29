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
    
    % Use shared helper for FFT-space corr2 in normalized Q-space
    [similarity_score, qpi_output, qpi_gt] = fftspace_corr2(output_kernel, gt_kernel);
    
    % Alternative metrics could be:
    % MSE = mean((qpi_output(:) - qpi_gt(:)).^2);
    % SSIM = ssim(qpi_output, qpi_gt);
end 