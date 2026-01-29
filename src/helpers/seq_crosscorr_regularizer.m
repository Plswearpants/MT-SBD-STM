function [ r ] = seq_crosscorr_regularizer(X_all_flat, gamma, sigma)
    %SEQ_CROSSCORR_REGULARIZER Creates a sequential cross-correlation regularizer
    %   r = seq_crosscorr_regularizer(X_all_flat, gamma, sigma)
    %
    %   Inputs:
    %       X_all_flat: flattened sum of all activation maps (required)
    %       gamma:      regularization parameter (default: 1e-3)
    %       sigma:      Gaussian smoothing parameter (default: 5.0)
    %
    %   Outputs:
    %       r:          struct containing cost, gradient, and Hessian functions
    %
    %   The regularizer penalizes high correlation between the current activation
    %   map (X) and the sum of all other activation maps (X_all_flat - X).

    if nargin < 1 || isempty(X_all_flat)
        error('X_all_flat must be provided as the flattened sum of all activation maps');
    end
    
    if nargin < 2 || isempty(gamma)
        gamma = 1e-3;
    end
    if nargin < 3 || isempty(sigma)
        sigma = 3.0;  % Default sigma value
    end

    % Validate input types and sizes
    if ~isnumeric(X_all_flat) || ~ismatrix(X_all_flat)
        error('X_all_flat must be a 2D matrix');
    end
    if ~isnumeric(gamma) || ~isscalar(gamma) || gamma < 0
        error('gamma must be a non-negative scalar');
    end
    if ~isnumeric(sigma) || ~isscalar(sigma) || sigma <= 0
        error('sigma must be a positive scalar');
    end

    % Store parameters in struct
    r.gamma = gamma;
    r.sigma = sigma;
    r.X_all_flat = X_all_flat;
    
    % Store functions in struct
    r.cost = @(X, store) cost(X, r, store);
    r.grad = @(X, store) grad(X, r, store);
    r.hess = @(X, eta, store) hess(X, eta, r, store);
end

function [corr_val, x1_aligned, x2_aligned] = simple_crosscorr(x1, x2, sigma)
    % Simple cross-correlation with alignment
    % Inputs:
    %   x1, x2: input 2D maps
    %   sigma: Gaussian smoothing parameter
    % Outputs:
    %   corr_val: maximum correlation value
    %   x1_aligned, x2_aligned: aligned maps
    
    % Center the maps
    x1 = x1 - mean(x1(:));
    x2 = x2 - mean(x2(:));
    
    % Simple Gaussian smoothing
    window_size = ceil(2 * sigma);
    [x, y] = meshgrid(-window_size:window_size);
    gaussian_kernel = exp(-(x.^2 + y.^2)/(2*sigma^2));
    gaussian_kernel = gaussian_kernel / sum(gaussian_kernel(:));
    
    % Apply smoothing
    x1_smooth = conv2(x1, gaussian_kernel, 'same');
    x2_smooth = conv2(x2, gaussian_kernel, 'same');
    
    % Compute cross-correlation
    c = normxcorr2(x1_smooth, x2_smooth);
    
    % Find peak correlation
    [max_corr, max_idx] = max(c(:));
    [y_peak, x_peak] = ind2sub(size(c), max_idx);
    
    % Calculate offset
    y_offset = y_peak - size(x1, 1);
    x_offset = x_peak - size(x1, 2);
    
    % Align maps
    x2_aligned = circshift(x2, [-y_offset, -x_offset]);
    x1_aligned = x1;  % x1 is reference
    
    corr_val = max_corr;
end

function [ rx, store ] = cost(X, r, store)
    % Compute the cross-correlation cost between the current activation X
    % and the sum of all other activations (X_all_flat - X)
    
    % Initialize store if not provided
    if nargin < 4 || isempty(store)
        store = struct();
    end
    
    % Get the current activation map
    x_current = X;
    norm_current = norm(x_current(:));
    
    % Skip if current activation is zero
    if norm_current < eps
        rx = 0;
        return;
    end
    
    % Check if we have cached values
    if ~isfield(store, 'corr_val') || ~isfield(store, 'x_current_aligned') || ~isfield(store, 'x_other_aligned')
        % Compute correlation with flattened sum of all activations
        X_all_flat = r.X_all_flat;
        X_other = X_all_flat - X;  % Sum of all other activations
        [corr_val, x_current_aligned, x_other_aligned] = simple_crosscorr(x_current, X_other, r.sigma);
        
        % Store computed values
        store.corr_val = corr_val;
        store.x_current_aligned = x_current_aligned;
        store.x_other_aligned = x_other_aligned;
        store.norm_x = norm(x_current_aligned(:));
    end
    
    % Compute cost using cached values
    rx = r.gamma * store.corr_val^2; % square to make it more sensitive to larger correlations
end

function [ rg, store ] = grad(X, r, store)
    % Compute gradient of the cross-correlation regularizer
    % The gradient is computed with respect to the current activation X
    
    % Initialize store if not provided
    if nargin < 4 || isempty(store)
        store = struct();
    end
    
    % Get the current activation map
    x_current = X;
    norm_current = norm(x_current(:));
    
    % Skip if current activation is zero
    if norm_current < eps
        rg = zeros(size(X,1), size(X,2));
        return;
    end
    
    % Check if we have cached values
    if ~isfield(store, 'corr_val') || ~isfield(store, 'x_current_aligned') || ~isfield(store, 'x_other_aligned')
        % If not, compute them
        [~, store] = cost(X, r, store);
    end
    
    % Compute gradient using cached values
    rg = 2*r.gamma*store.corr_val * ...
        (store.x_other_aligned/store.norm_x - ...
        store.corr_val*store.x_current_aligned/store.norm_x^2);
end

function [ rh, store ] = hess(X, eta, r, store)
    % Compute Hessian-vector product of the cross-correlation regularizer
    % eta is the direction vector for the current activation X
    
    % Initialize store if not provided
    if nargin < 5 || isempty(store)
        store = struct();
    end
    
    % Get the current activation map
    x_current = X;
    norm_current = norm(x_current(:));
    
    % Skip if current activation is zero
    if norm_current < eps
        rh = zeros(size(eta));
        return;
    end
    
    % Check if we have cached values
    if ~isfield(store, 'corr_val') || ~isfield(store, 'x_current_aligned') || ~isfield(store, 'x_other_aligned')
        % If not, compute them
        [~, store] = cost(X, r, store);
    end
    
    % Compute Hessian-vector product using cached values
    norm_x_sq = store.norm_x^2;
    
    % First term: 2*gamma*corr_val * (x_other_aligned/norm_x)
    term1 = 2*r.gamma*store.corr_val * store.x_other_aligned/store.norm_x;
    
    % Second term: -2*gamma*corr_val^2 * (x_current_aligned/norm_x_sq)
    term2 = -2*r.gamma*store.corr_val^2 * store.x_current_aligned/norm_x_sq;
    
    % Third term: 2*gamma * (x_other_aligned'*eta/norm_x - corr_val*x_current_aligned'*eta/norm_x_sq) * x_current_aligned
    inner_prod1 = sum(sum(store.x_other_aligned.*eta));
    inner_prod2 = sum(sum(store.x_current_aligned.*eta));
    term3 = 2*r.gamma * (inner_prod1/store.norm_x - store.corr_val*inner_prod2/norm_x_sq) * store.x_current_aligned;
    
    % Combine terms
    rh = term1 + term2 + term3;
end 