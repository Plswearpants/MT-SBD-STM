function [ Xsol, info ] = Xsolve_FISTA_multikernel( Y, A, lambda, mu, varargin )
    %XSOLVE_FISTA_MULTIKERNEL   Solve for X using FISTA method with multiple kernels
    % Note that this function assumes only one energy layer.
    %   - Core usage:
    %       [ Xsol, info ] = Xsolve_FISTA_multikernel( Y, A, lambda, mu )
    %
    %   - Optional variables:
    %       [ ... ] = Xsolve_FISTA_multikernel( ... , Xinit, Xpos, getbias )
    %       Xinit:      initial value for X
    %       Xpos:       constrain X to be a positive solution
    %       getbias:    extract constant bias as well as X
    %
    
    % Initialize variables and function handles:
    fpath = fileparts(mfilename('fullpath'));
    load([fpath '/../../config/Xsolve_config.mat']); %#ok<*LOAD>
    g = huber(mu);
    num_kernels = size(A,3);
    m = size(Y);
    % define size of X as spatial size + number of kernels
    X_size = [m, num_kernels];

    %% Checking arguments:
    nvararg = numel(varargin);
    if nvararg > 3
        error('Too many input arguments.');
    end

    idx = 1; X = zeros(X_size); b = 0; % note that bias is now a scalar
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
    
    %% Iterate:    
    t=1; W = X; u = b;
    costs = NaN(MAXIT,2);
    doagain = true;  it = 0;  count = 0;
    while doagain
        it = it + 1;
        % Gradients and Hessians:
        grad_fW = zeros(X_size); grad_fu = 0;
        Ri = zeros(m);
        for k = 1:num_kernels
            Ri = Ri + convfft2(A(:,:,k), W(:,:,k));
        end
        Ri = Ri + u - Y;

        for k = 1:num_kernels
            grad_fW(:,:,k) = grad_fW(:,:,k) + convfft2( A(:,:,k), Ri, 1 );
        end
        grad_fu = sum(Ri(:));
        
        %???
        A_fft = fft2(A, m(1), m(2));
        R_A = sum(abs(A_fft).^2, 3);


        % FISTA update
        L = max(R_A(:));
        X_ = g.prox(W - 1/L*grad_fW, lambda/L, xpos);
        t_ = (1+sqrt(1+4*t^2))/2;
        W = X_ + (t-1)/t_*(X_-X);
        if getbias
            b_ = u - grad_fu/(2*prod(m));  
            u = b_ + (t-1)/t_*(b_-b);
            b = b_;
        end
        X = X_; t = t_;

        % Check conditions to repeat iteration:
        p = zeros(m);
        for k = 1:num_kernels
            p = p + convfft2(A(:,:,k), X(:,:,k));
        end
        f = norm(p + b - Y, 'fro')^2/2;
        costs(it,1) = f;
        costs(it,2) = g.cost(X, lambda);

        tmp = grad_fW + grad_fu;
        delta = g.diffsubg(X, -tmp, lambda, xpos);
        delta = norm(delta(:))/sqrt(prod(X_size));
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
    Xsol.f = sum(costs(it,:));
    info.numit = it;
    info.delta = delta;
    info.costs = costs(1:it,:);
end
