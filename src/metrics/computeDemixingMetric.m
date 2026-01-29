function [demixing_score, overlap_matrix] = computeDemixingMetric(X, metric_type, binary_mode)
    % Compute separation metric for activation maps using IoU or Dice
    %
    % Input:
    %   X: activation maps [height x width x num_kernels]
    %   metric_type: string specifying the metric to use:
    %       'iou'  - Intersection over Union (default)
    %       'dice' - Dice coefficient
    %   binary_mode: boolean, whether to binarize activations (default: false)
    %
    % Output:
    %   demixing_score: scalar value representing separation quality
    %                    (1 = perfect separation, 0 = complete overlap)
    %   overlap_matrix: matrix of overlap measures between activations
    
    if nargin < 2
        metric_type = 'iou';
    end
    if nargin < 3
        binary_mode = false;
    end
    
    num_kernels = size(X, 3);
    overlap_matrix = zeros(num_kernels, num_kernels);
    
    % Normalize activation maps to [0,1] if not binary mode
    if ~binary_mode
        % Normalize each activation map independently
        for i = 1:num_kernels
            X_slice = X(:,:,i);
            X(:,:,i) = (X_slice - min(X_slice(:))) / (max(X_slice(:)) - min(X_slice(:)) + eps);
        end
    else
        X = X > 0;
    end
    
    % Compute overlap between all pairs
    for i = 1:num_kernels
        X1 = X(:,:,i);
        for j = i+1:num_kernels
            X2 = X(:,:,j);
            
            switch lower(metric_type)
                case 'iou'
                    % Intersection over Union
                    intersection = sum(min(X1(:), X2(:)));
                    union = sum(max(X1(:), X2(:)));
                    overlap = intersection / (union + eps);
                    
                case 'dice'
                    % Dice coefficient
                    intersection = sum(X1(:) .* X2(:));
                    overlap = 2 * intersection / (sum(X1(:)) + sum(X2(:)) + eps);
                    
                otherwise
                    error('Unknown metric type: use ''iou'' or ''dice''');
            end
            
            % Store in symmetric matrix
            overlap_matrix(i,j) = overlap;
            overlap_matrix(j,i) = overlap;
        end
    end
    
    % Calculate demixing score as 1 - average overlap
    if num_kernels > 1
        avg_overlap = sum(overlap_matrix(:)) / (num_kernels * (num_kernels - 1));
    else
        avg_overlap = 0;  % Single kernel case
    end
    demixing_score = 1 - avg_overlap;
end 