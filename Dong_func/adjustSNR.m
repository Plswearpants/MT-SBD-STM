function Y_new = adjustSNR(Y, A0, X0, current_SNR, target_SNR)
    % Adjust noise level in Y to achieve target SNR while preserving A0, X0
    %
    % Inputs:
    %   Y: Current observation with noise
    %   A0, X0: Original kernel and activation
    %   current_SNR: Current SNR of Y
    %   target_SNR: Desired SNR
    %
    % Output:
    %   Y_new: Observation with adjusted noise level
    
    % Generate clean signal
    Y_clean = zeros(size(Y));
    for k = 1:size(X0,3)
        Y_clean = Y_clean + convfft2(A0{k}, X0(:,:,k));
    end
    
    % Extract current noise
    current_noise = Y - Y_clean;
    
    % Calculate scaling factor for noise
    % SNR = var(signal)/var(noise)
    % target_SNR = var(signal)/(alpha^2 * var(current_noise))
    % Therefore: alpha = sqrt(current_SNR/target_SNR)
    alpha = sqrt(current_SNR/target_SNR);
    
    % Scale noise and add back to clean signal
    Y_new = Y_clean + alpha * current_noise;
end 