function quality_factor = compute_kernel_quality_factors(ground_truth_kernel, input_kernel, noise_level)
    % Compute quality factor between ground truth and input kernels with noise normalization
    %
    % Inputs:
    %   ground_truth_kernel: reference kernel
    %   input_kernel: kernel to be evaluated
    %   noise_level: (optional) noise level for normalization, default is estimated from data
    %
    % Output:
    %   quality_factor: normalized quality metric (lower is better)
    
    % Ensure both kernels are column vectors
    ground_truth_kernel = ground_truth_kernel(:);
    input_kernel = input_kernel(:);
    
    % Estimate noise level if not provided
    if nargin < 3 || isempty(noise_level)
        % Estimate noise from ground truth kernel
        % Using median absolute deviation (MAD) for robust estimation
        diff_kernel = diff(ground_truth_kernel);
        noise_level = 1.4826 * median(abs(diff_kernel - median(diff_kernel)));
        fprintf('Estimated noise level: %.2e\n', noise_level);
    end
    
    % Add small epsilon to prevent division by zero
    eps_noise = max(noise_level, 1e-12);
    
    % Compute normalized quality factor
    diff_variance = var(ground_truth_kernel - input_kernel);
    quality_factor = diff_variance / (eps_noise^2);
    
    % Print diagnostics
    fprintf('Quality Factor Components:\n');
    fprintf('  Difference Variance: %.2e\n', diff_variance);
    fprintf('  Noise Level^2: %.2e\n', eps_noise^2);
    fprintf('  Normalized Quality Factor: %.2e\n', quality_factor);
end

