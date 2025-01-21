function visualizeResults_old(Y, A0, Aout, X0, Xout, bout, extras)
    % Visualize results from SBD-STM reconstruction
    %
    % Inputs:
    %   Y: Original observation
    %   A0: Ground truth kernels (cell array)
    %   Aout: Reconstructed kernels (cell array)
    %   X0: Ground truth activation maps
    %   Xout: Reconstructed activation maps
    %   bout: Bias terms [num_kernels x 1]
    %   extras: Struct containing convergence metrics
    
    % Convert cell arrays to matrices
    extras.phase1.activation_metrics = old2new_metric(extras.phase1.activation_metrics);
    extras.phase2.activation_metrics = old2new_metric(extras.phase2.activation_metrics);
    
    % call the function to visualize the results
    visualizeResults(Y, A0, Aout, X0, Xout, bout, extras);
end 
