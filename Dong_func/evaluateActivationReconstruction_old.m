function [metrics, aligned_maps] = evaluateActivationReconstruction_old(X0, Xout, kernel_size, visualize)
    % Evaluate reconstruction quality by first aligning activation maps and then
    % computing similarity with density-adaptive Gaussian filtering
    %
    % Inputs:
    %   X0: original activation maps [height x width x num_kernels]
    %   Xout: reconstructed activation maps [height x width x num_kernels]
    %   kernel_size: size of each kernel [num_kernels x 2]
    %   visualize: (optional) boolean for visualization (must be 1 or 0)
    %
    % Outputs:
    %   metrics: struct containing alignment and similarity metrics per kernel
    %   aligned_maps: aligned activation maps and filtered results
    
    if nargin < 4
        visualize = false;
    end
    
    % Step 1: Align activation maps
    [Xout_aligned, offsets, align_quality] = alignActivationMaps(X0, Xout, kernel_size);
    
    % Step 2: Compute similarity with density-adaptive Gaussian
    [similarities, filtered_maps] = computeActivationSimilarity(X0, Xout_aligned, kernel_size, visualize);
    
    % Combine metrics for each kernel
    num_kernels = size(X0, 3);
    for k = 1:num_kernels
        metrics(k).alignment = align_quality(k);     % Alignment quality metrics
        metrics(k).offset = offsets(k,:);            % Applied offset
        metrics(k).similarity = similarities(k);      % Final similarity score
        metrics(k).density = filtered_maps(k).density;% Activation density
        metrics(k).L = filtered_maps(k).L;           % Characteristic length
        metrics(k).sigma = filtered_maps(k).sigma;   % Used Gaussian width
    end
    
    % Store aligned and filtered maps
    aligned_maps.Xout_aligned = Xout_aligned;
    aligned_maps.filtered = filtered_maps;
    
    % Print comprehensive results
    fprintf('\nFinal Evaluation Results:\n');
    fprintf('========================\n');
    for k = 1:num_kernels
        fprintf('\nKernel %d:\n', k);
        fprintf('Alignment:\n');
        fprintf('  Pattern Match: %.3f  (1 = perfect match)\n', metrics(k).alignment.primary_peak);
        fprintf('  Uniqueness: %.3f     (1 = single solution)\n', 1 - metrics(k).alignment.peak_ratio);
        fprintf('  Precision: %.3f      (>1 = precise alignment)\n', metrics(k).alignment.peak_sharpness);
        fprintf('  Offset: [%d, %d]\n', metrics(k).offset(1), metrics(k).offset(2));
        fprintf('Similarity:\n');
        fprintf('  Score: %.3f          (1 = perfect match)\n', metrics(k).similarity);
        fprintf('  Density: %.2e\n', metrics(k).density);
        fprintf('  L: %.1f pixels       (isolation length)\n', metrics(k).L);
        fprintf('  σ: %.1f pixels       (filter width)\n', metrics(k).sigma);
    end
end
