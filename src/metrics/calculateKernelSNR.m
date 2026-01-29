function SNR = calculateKernelSNR(A0, A0_noiseless)
    % Input validation
    validateattributes(A0, {'cell'}, {'vector'});
    validateattributes(A0_noiseless, {'cell'}, {'vector'});
    assert(length(A0) == length(A0_noiseless), 'A0 and A0_noiseless must have same number of kernels');
    
    num_kernels = length(A0);
    avg_var = 0;
    avg_noise_var = 0;
    
    % Calculate average variance of clean signal and noise
    for n = 1:num_kernels
        avg_var = avg_var + var(A0_noiseless{n}, 0, "all");
        noise = A0{n} - A0_noiseless{n};
        avg_noise_var = avg_noise_var + var(noise, 0, "all");
    end
    
    % Average across kernels
    avg_var = avg_var / num_kernels;
    avg_noise_var = avg_noise_var / num_kernels;
    
    % Calculate SNR as ratio of signal variance to noise variance
    SNR = avg_var / avg_noise_var;
end