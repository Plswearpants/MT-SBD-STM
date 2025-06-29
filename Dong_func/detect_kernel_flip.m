function [is_flipped, quality_scores, Xout, Aout] = detect_kernel_flip(X0, Xout, A0, Aout, skip_flip_check)
    % Detect if kernels are flipped by comparing quality metrics, if flipped, flip the kernels and activations, so the output is in the order that mathces the ground truth.
    %
    % Inputs:
    %   X0: Ground truth activations {1x1 cell} or direct [height x width x 2]
    %   Xout: Output activations {1x1 cell} or direct [height x width x 2]
    %   A0: Ground truth kernels {1x1 cell} or direct {1x2 cell}
    %   Aout: Output kernels {1x1 cell} or direct {1x2 cell}
    %   skip_flip_check: (Optional) Boolean to skip flip quality calculations
    %
    % Outputs:
    %   is_flipped: Boolean indicating if kernels are flipped
    %   quality_scores: Structure containing quality scores for both cases
    
    % Set default value for skip_flip_check if not provided
    if nargin < 5
        skip_flip_check = false;
    end
    
    % Handle cell array inputs
    if iscell(X0) && ~isempty(X0), X0 = X0{1}; end
    if iscell(Xout) && ~isempty(Xout), Xout = Xout{1}; end
    if iscell(A0) && ~isempty(A0), A0 = A0{1}; end
    if iscell(Aout) && ~isempty(Aout), Aout = Aout{1}; end
    
    % Get kernel size from A0
    if ~isempty(A0) && iscell(A0) && numel(A0) >= 2
        kernel_size = [size(A0{1}); size(A0{2})];
    else
        error('Invalid kernel structure');
    end

    % Case 1: No flip - original order
    [activation_similarity, kernel_similarity] = computeQualityMetrics(X0, Xout, A0, Aout, kernel_size);
    
    % Store quality scores
    quality_scores = struct();
    quality_scores.no_flip = struct('activation_similarity', activation_similarity, ...
                                  'kernel_similarity', kernel_similarity);
    
    if skip_flip_check
        is_flipped = false;
    else
        % Case 2: Flipped case - swap kernels and activations
        Xout_flip = cat(3, Xout(:,:,2), Xout(:,:,1));
        Aout_flip = {Aout{2}, Aout{1}};  % Swap kernels
        
        [activation_similarity_flipped, kernel_similarity_flipped] = computeQualityMetrics(X0, Xout_flip, A0, Aout_flip, kernel_size);
        
        % Add flipped quality scores
        quality_scores.flipped = struct('activation_similarity', activation_similarity_flipped, ...
                                      'kernel_similarity', kernel_similarity_flipped);

        % Determine if kernels should be flipped - flipped if all entries in quality are greater in flipped case
        is_flipped = all(activation_similarity_flipped > activation_similarity) && all(kernel_similarity_flipped > kernel_similarity) && all(activation_similarity_flipped > 0.8) && all(kernel_similarity_flipped > 0.8);
        
        % Only print information if kernels are flipped
        if is_flipped
            fprintf('Kernels are FLIPPED (activation similarity improvement: %.3f -> %.3f , kernel similarity improvement: %.3f -> %.3f)\n', ...
                quality_scores.no_flip.activation_similarity, quality_scores.flipped.activation_similarity, ...
                quality_scores.no_flip.kernel_similarity, quality_scores.flipped.kernel_similarity);
            % flip the kernels and activations
            Xout = Xout_flip;
            Aout = Aout_flip;
        end
    end
end
