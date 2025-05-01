function [ Xsol, info ] = Xsolve_FISTA_ALL( Y, A, lambda, mu, varargin )
%XSOLVE_FISTA   Solve for X using FISTA method
%   - Core usage:
%       [ Xsol, info ] = Xsolve_FISTA( Y, A, lambda, mu )
%
%   - Optional variables:
%       [ ... ] = Xsolve_FISTA( ... , Xinit, Xpos, getbias )
%       Xinit:      initial value for X
%       Xpos:       constrain X to be a positive solution
%       getbias:    extract constant bias as well as X
%

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
    
    % Get number of types from A
    A_size = size(A);
    if numel(A_size) > 3
        t = A_size(1); % Number of types
    else
        t = 1; % Single type
        A = reshape(A, [t, A_size]); % Reshape to new dimension format
    end
  
    %% Checking arguments:
    nvararg = numel(varargin);
    if nvararg > 3
        error('Too many input arguments.');
    end

    idx = 1; 
    if nvararg >= idx && ~isempty(varargin{idx})
        X = varargin{idx}.X;
        b = varargin{idx}.b;
    else
        X = zeros([t, m]); % Initialize X with type dimension
        b = zeros(n,1);
    end

    idx = 2; xpos = false;
    if nvararg >= idx && ~isempty(varargin{idx})
        xpos = varargin{idx};
    end
    
    idx = 3; getbias = false;
    if nvararg >= idx && ~isempty(varargin{idx})
        getbias = varargin{idx};
    end

    
    
    %% Iterate:    
    t_iter = 1; W = X; u = b;
    costs = NaN(MAXIT,2);
    doagain = true;  it = 0;  count = 0;
    while doagain
        it = it + 1;
        % Gradients and Hessians:
        grad_fW = zeros(size(X)); grad_fu = zeros(n,1); R_A = zeros(m);
        for i = 1:n     % sum up over slices
            Ri = zeros(m);
            for j = 1:t  % sum up over types
                Ri = Ri + convfft2(A(j,:,:,i), W(j,:,:));
            end
            Ri = Ri + u(i) - Y(:,:,i);
            
            for j = 1:t  % compute gradients for each type
                grad_fW(j,:,:) = grad_fW(j,:,:) + convfft2(A(j,:,:,i), Ri, 1);
            end
            
            grad_fu(i) = sum(Ri(:));
            
            % Compute R_A (used for step size)
            for j = 1:t
                R_A = R_A + abs(fft2(A(j,:,:,i),m(1),m(2))).^2;
            end
        end

        % FISTA update
        L = max(R_A(:));
        
        % Update X for each type
        X_ = zeros(size(X));
        for j = 1:t
            X_(j,:,:) = g.prox(W(j,:,:) - 1/L*grad_fW(j,:,:), lambda/L, xpos);
        end
        
        t_ = (1+sqrt(1+4*t_iter^2))/2;
        
        % Momentum step for each type
        for j = 1:t
            W(j,:,:) = X_(j,:,:) + (t_iter-1)/t_*(X_(j,:,:)-X(j,:,:));
        end
        
        if getbias
            b_ = u - grad_fu/(2*prod(m)*sqrt(n));  
            u = b_ + (t_iter-1)/t_*(b_-b);
            b = b_;
        end
        X = X_; t_iter = t_;

        % Check conditions to repeat iteration:
        f = 0;
        for i = 1:n
            total_conv = zeros(m);
            for j = 1:t
                total_conv = total_conv + convfft2(A(j,:,:,i), reshape(X(j,:,:), m));
            end
            f = f + norm(total_conv + b(i) - Y(:,:,i), 'fro')^2/2;
        end
        costs(it,1) = f;
        
        % Sum sparsity costs over all types
        reg_cost = 0;
        for j = 1:t
            reg_cost = reg_cost + g.cost(X(j,:,:), lambda);
        end
        costs(it,2) = reg_cost;

        % Compute convergence criterion
        delta = 0;
        for j = 1:t
            tmp = grad_fW(j,:,:);
            tmp_delta = g.diffsubg(X(j,:,:), -tmp, lambda, xpos);
            delta = max(delta, norm(tmp_delta(:))/sqrt(prod(m)));
        end
        
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
    Xsol.W = W;         % For compatibility with pdNCG
    Xsol.f = sum(costs(it,:));
    info.numit = it;
    info.delta = delta;
    info.costs = costs(1:it,:);
end
