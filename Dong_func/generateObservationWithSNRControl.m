function [Y, Y_clean, A0, num_kernels, kernel_size] = generateObservationWithSNRControl(X0, A0_noiseless, SNR)
    % Generate observation and kernels with controlled SNR
    % Inputs:
    %   X0: Activation maps (3D array)
    %   A0_noiseless: Cell array of noiseless kernels
    %   SNR: Target signal-to-noise ratio
    % Outputs:
    %   Y: Noisy observation
    %   Y_clean: Clean observation before noise addition
    %   A0: Cell array of noisy kernels with same SNR as observation
    
    % Generate clean observation
    Y_clean = zeros(size(X0, 1:2));
    for k = 1:size(X0, 3)
        Y_clean = Y_clean + convfft2(A0_noiseless{k}, X0(:,:,k));
    end
    
    % Add noise to match target SNR for observation
    eta = var(Y_clean, 0, "all") / SNR;
    Y = Y_clean + sqrt(eta) * randn(size(Y_clean));
    
    % Calculate average variance of kernels
    num_kernels = length(A0_noiseless);
    kernel_size = zeros(num_kernels,2);
    total_var = 0;
    for n = 1:num_kernels
        A0_noiseless{n} = proj2oblique(A0_noiseless{n});
        kernel_size(n,:) = size(A0_noiseless{n});
        total_var = total_var + var(A0_noiseless{n}, 0, "all");
    end
    avg_var = total_var / num_kernels;
    
    % Apply universal noise to kernels based on average variance
    eta_kernel = avg_var / SNR;  % Universal noise variance
    A0 = cell(1, num_kernels);
    for n = 1:num_kernels
        A0{n} = A0_noiseless{n} + sqrt(eta_kernel) * randn(size(A0_noiseless{n}));
        A0{n} = proj2oblique(A0{n});  % Project back to sphere after adding noise
    end
end 