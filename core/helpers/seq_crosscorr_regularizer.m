function [ r ] = seq_crosscorr_regularizer( gamma )
    if nargin < 1 || isempty(gamma)
        gamma = 1e-3;
    end

    r.cost = @(X, current_idx) cost(X, current_idx, gamma);
    r.grad = @(X, current_idx) grad(X, current_idx, gamma);
end

function [ cost ] = cost(X, current_idx, gamma)
    % X is a 3D matrix where each slice X(:,:,i) represents an activation map
    % current_idx is the index of the activation being updated
    n = size(X, 3);
    cost = 0;
    
    % Get the current activation map
    x_current = X(:,:,current_idx);
    x_current = x_current - mean(x_current(:));
    norm_current = norm(x_current(:));
    
    % Skip if current activation is zero
    if norm_current < eps
        return;
    end
    
    % Create a 3D array with current activation repeated
    X_current = repmat(x_current, [1, 1, n-1]);
    X_others = X(:,:,setdiff(1:n, current_idx));
    
    % Compute kernel size based on activation size
    [m, n] = size(x_current);
    kernel_size = repmat([m, n], [n-1, 1]);
    
    % Compute similarity using the same function as evaluation
    [similarities, ~] = computeActivationSimilarity(X_current, X_others, kernel_size);
    
    % Sum up similarities (we want to minimize correlation)
    cost = gamma * sum(similarities.^2);
end

function [ grad ] = grad(X, current_idx, gamma)
    % Compute gradient of the cross-correlation regularizer
    % Only compute gradient for the current activation
    n = size(X, 3);
    grad = zeros(size(X,1), size(X,2));
    
    % Get the current activation map
    x_current = X(:,:,current_idx);
    x_current = x_current - mean(x_current(:));
    norm_current = norm(x_current(:));
    
    % Skip if current activation is zero
    if norm_current < eps
        return;
    end
    
    % Create a 3D array with current activation repeated
    X_current = repmat(x_current, [1, 1, n-1]);
    X_others = X(:,:,setdiff(1:n, current_idx));
    
    % Compute kernel size based on activation size
    [m, n] = size(x_current);
    kernel_size = repmat([m, n], [n-1, 1]);
    
    % Compute similarity and get filtered maps
    [similarities, filtered_maps] = computeActivationSimilarity(X_current, X_others, kernel_size);
    
    % Compute gradient using filtered maps
    for k = 1:length(similarities)
        % Get filtered maps for this pair
        X_current_filtered = filtered_maps(k).X0;
        X_other_filtered = filtered_maps(k).Xout;
        
        % Compute gradient component
        grad = grad + 2*gamma*similarities(k) * ...
            (X_other_filtered/norm(X_current_filtered(:)) - ...
            similarities(k)*X_current_filtered/norm(X_current_filtered(:))^2);
    end
    
    % Apply the same Gaussian filter used in computeActivationSimilarity
    sigma = filtered_maps(1).sigma;
    window_size = ceil(3 * sigma);
    [x, y] = meshgrid(-window_size:window_size);
    gaussian_kernel = exp(-(x.^2 + y.^2)/(2*sigma^2));
    gaussian_kernel = gaussian_kernel / sum(gaussian_kernel(:));
    
    % Smooth the gradient
    grad = conv2(grad, gaussian_kernel, 'same');
end 