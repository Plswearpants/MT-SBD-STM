function [ Xsol, info ] = Xsolve_FISTA_test( Y, A, lambda, mu, X_all_flat, gamma, varargin )
%XSOLVE_FISTA_TEST   Test version of Xsolve_FISTA with sequential cross-correlation regularizer
%   - Core usage:
%       [ Xsol, info ] = Xsolve_FISTA_test( Y, A, lambda, mu, X_all_flat, gamma )
%
%   - Optional variables:
%       [ ... ] = Xsolve_FISTA_test( ... , Xinit, Xpos, getbias )
%       Xinit:      initial value for X
%       Xpos:       constrain X to be a positive solution
%       getbias:    extract constant bias as well as X

    % Validate X_all_flat
    if isempty(X_all_flat) || ~isnumeric(X_all_flat) || ~ismatrix(X_all_flat)
        error('X_all_flat must be a non-empty 2D matrix');
    end

    % Validate gamma
    if ~isnumeric(gamma) || ~isscalar(gamma) || gamma < 0
        error('gamma must be a non-negative scalar');
    end

    % Initialize variables and function handles:
    fpath = fileparts(mfilename('fullpath'));
    addpath([fpath '/helpers']);
    load([fpath '/../config/Xsolve_config.mat']); %#ok<*LOAD>
    g = huber(mu);

    m = size(Y);
    if (numel(m) > 2)
        n = m(3); m = m(1:2);
    else
        n = 1;
    end

    %% Checking arguments:
    nvararg = numel(varargin);
    if nvararg > 3
        error('Too many input arguments.');
    end

    idx = 1; X = zeros(m); b = zeros(n,1);
    if nvararg >= idx && ~isempty(varargin{idx})
        X = varargin{idx}.X;
        b = varargin{idx}.b;
    end

    idx = 2; xpos = false;
    if nvararg >= idx && ~isempty(varargin{idx})
        xpos = varargin{idx};
    end
    
    idx = 3; getbias = false;
    if nvararg >= idx && ~isempty(varargin{idx})
        getbias = varargin{idx};
    end

    % Initialize regularizer with X_all_flat and gamma
    r = seq_crosscorr_regularizer(X_all_flat, gamma);

    % Initialize W and u of the current kernel
    W = X;
    u = b;
    t = 1;
    costs = zeros(MAXIT, 3);  % [data_fidelity, regularization, cross_correlation]
    delta = inf;
    count = 0;
    it = 0;
    doagain = true;

    % Initialize store for cross-correlation regularizer
    store = struct();

    while doagain
        it = it + 1;
        grad_fW = zeros(m);
        grad_fu = zeros(n,1);
        R_A = zeros(m);

        % Compute gradients for each slice
        for i = 1:n
            Ri = convfft2(A(:,:,i), W) + u(i) - Y(:,:,i);
            grad_fW = grad_fW + convfft2( A(:,:,i), Ri, 1 );
            grad_fu(i) = sum(Ri(:));
            R_A = R_A + abs(fft2(A(:,:,i),m(1),m(2))).^2;
        end

        % Add scaled cross-correlation gradient
        if gamma > 0
            [cross_corr_grad, store] = r.grad(W, store);
            % Ensure gradient is finite and scale by gamma
            if any(~isfinite(cross_corr_grad(:)))
                warning('Non-finite values in cross-correlation gradient, setting to zero');
                cross_corr_grad(~isfinite(cross_corr_grad)) = 0;
            end
            grad_fW = grad_fW + cross_corr_grad;  
        end

        % FISTA update
        L = max(R_A(:));
        X_ = g.prox(W - 1/L*grad_fW, lambda/L, xpos);
        % Ensure X_ is finite
        if any(~isfinite(X_(:)))
            warning('Non-finite values in X_, setting to zero');
            X_(~isfinite(X_)) = 0;
        end
        t_ = (1+sqrt(1+4*t^2))/2;
        W = X_ + (t-1)/t_*(X_-X);
        if getbias
            b_ = u - grad_fu/(2*prod(m)*sqrt(n));  
            % Ensure b_ is finite
            if any(~isfinite(b_(:)))
                warning('Non-finite values in b_, setting to zero');
                b_(~isfinite(b_)) = 0;
            end
            u = b_ + (t-1)/t_*(b_-b);
            b = b_;
        end
        X = X_; t = t_;

        % Compute costs
        f = 0;
        for i = 1:n
            f = f + norm(convfft2(A(:,:,i), reshape(X, m)) + b(i) - Y(:,:,i), 'fro')^2/2;
        end
        costs(it,1) = f;
        costs(it,2) = g.cost(X, lambda);  % regularization

        if gamma > 0
            [cross_corr_cost, store] = r.cost(X, store);
            % Ensure cost is finite
            if ~isfinite(cross_corr_cost)
                warning('Non-finite cross-correlation cost, setting to zero');
                cross_corr_cost = 0;
            end
            costs(it,3) = cross_corr_cost;  % cross-correlation
        else
            costs(it,3) = 0;
        end

        % Check conditions to repeat iteration:
        tmp = grad_fW;
        for i = 1:n
            tmp = tmp + grad_fu(i);
        end
        delta = g.diffsubg(X, -tmp, lambda, xpos);
        delta = norm(delta(:))/sqrt(prod(m));
        if delta < EPSILON
            count = count+1;
        else
            count = 0;
        end
        doagain = count < 10 && (it < MAXIT);
    end
    
    % Return solution:
    Xsol.X = X;
    Xsol.b = b;
    Xsol.W = W;         % dummy variable for compatibility with pdNCG.
    Xsol.f = sum(costs(it,:));  % Total cost including all terms
    info.numit = it;
    info.delta = delta;
    info.costs = costs(1:it,:);  % Now includes cross-correlation costs
end 