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
    if iscell(A0) && ~isempty(A0) && numel(A0) <= 1, A0 = A0{1}; end
    if iscell(Aout) && ~isempty(Aout), Aout = Aout{1}; end
    
    % Get kernel size from A0 (generalized for any number of kernels)
    if ~isempty(A0) && iscell(A0)
        num_kernels = numel(A0);
        kernel_size = zeros(num_kernels, 2);
        for k = 1:num_kernels
            kernel_size(k,:) = size(A0{k});
        end
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
        % For more than 2 kernels, check all permutations
        nK = numel(A0);
        perms_to_check = perms(1:nK);
        best_score = -Inf;
        best_perm = 1:nK;
        for p = 1:size(perms_to_check,1)
            perm = perms_to_check(p,:);
            Xout_perm = Xout(:,:,perm);
            Aout_perm = Aout(perm);
            [as_perm, ks_perm] = computeQualityMetrics(X0, Xout_perm, A0, Aout_perm, kernel_size);
            score = mean(as_perm) + mean(ks_perm);
            if score > best_score
                best_score = score;
                best_perm = perm;
                best_as = as_perm;
                best_ks = ks_perm;
            end
        end
        % Save best permutation as 'flipped'
        quality_scores.flipped = struct('activation_similarity', best_as, ...
                                        'kernel_similarity', best_ks);
        is_flipped = ~isequal(best_perm, 1:nK);
        if is_flipped
            fprintf('Kernels are FLIPPED (best permutation: %s)\n', mat2str(best_perm));
            Xout = Xout(:,:,best_perm);
            Aout = Aout(best_perm);
        end
    end
end
