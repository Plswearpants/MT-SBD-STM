function [ A_phase_I, Aout, Xout, bout, extras ] = SBD_test( Y, k, params, dispfun, kernel_data )
%SBD Summary of this function goes here
%
%   PARAMS STRUCT:
%   ===============
%   The options struct should include the fields:
%       lambda1,  float > 0  : regularization parameter for Phase I
%       phase2,   bool       : whether to do Phase II (refinement) or not
%
%   IF phase2 == true, then the following fields should also be included:
%       kplus,    int > 0    : border padding (pixels) for sphere lifting
%       lambda2,  float > 0  : FINAL reg. param. value for Phase II
%
%       nrefine,  int >= 1   : number of refinements for Phase II.
%           Refinement 1 lifts the sphere and uses lambda1, successive
%           refinements decrease lambda down to lambda2;
%           i.e. if nrefine == 1, then no decrease in lambda is made.
%
%
%   Finally, two optional fields for the struct. These features are
%   automatically disabled if the fields are not included or are empty:
%
%       Xsolve,   string     : Pick which Xsolve to use--'FISTA' or
%       'pdNCG.'
%
%       xpos,     bool       :  Constrain X to have nonnegative entries
%           when running XSolve.
%
%       getbias,  bool       : Extract constant bias from observation.
%                                                                                                                      

%% Process input arguments
starttime = tic;
n = size(Y,3);

if nargin < 4 || isempty(dispfun)
    dispfun = @(Y,A,X,k,kplus,idx) 0;
end

lambda1 = params.lambda1;
if params.phase2
    kplus = params.kplus;
    lambda2 = params.lambda2;
    nrefine = params.nrefine;
end

if ~isfield(params, 'xpos') || isempty(params.xpos)
    xpos = false;
else
    xpos = params.xpos;
end


if ~isfield(params, 'getbias') || isempty(params.getbias)
    getbias = false;
else
    getbias = params.getbias;
end

if ~isfield(params, 'Xsolve') || isempty(params.Xsolve)
    Xsolve = 'FISTA';
else
    Xsolve = params.Xsolve;
end

%% PHASE I: First pass at BD
dispfun1 = @(A, X) dispfun(Y, A, X, k, [], 1);

fprintf('PHASE I: \n=========\n');
A = kernel_data;

[A, Xsol, info] = Asolve_Manopt( Y, A, lambda1, Xsolve, [], xpos, getbias, dispfun1);
extras.phase1.A = A;
extras.phase1.X = Xsol.X;
extras.phase1.b = Xsol.b;
extras.phase1.info = info;
A_phase_I = A ;
%% PHASE II: Lift the sphere and do lambda continuation
if params.phase2
    k2 = k + 2*kplus;
    dispfun2 = @(A, X) dispfun(Y, A, X, k2, 0, 1);

    A2 = zeros([k2 n]);
    A2(kplus(1)+(1:k(1)), kplus(2)+(1:k(2)), :) = A;
    X2sol = Xsol;
    %X2sol.X = circshift(Xsol.X,-kplus);
    %X2sol.W = circshift(Xsol.W,-kplus);
    % clear A Xsol;

    lambda = lambda1;
    score = zeros(2*kplus+1);
    fprintf('\n\nPHASE II: \n=========\n');
    lam2fac = (lambda2/lambda1)^(1/nrefine);
    i = 1;
    while i <= nrefine + 1
        fprintf('lambda = %.1e: \n', lambda);
        [A2, X2sol, info] = Asolve_Manopt( Y, A2, lambda, Xsolve, X2sol, xpos, getbias, dispfun2 );
        fprintf('\n');

        %Attempt to 'unshift" the a and x by taking the l1-norm over all k-contiguous elements:
        for tau1 = -kplus(1):kplus(1)
            ind1 = tau1+kplus(1)+1;
            for tau2 = -kplus(2):kplus(2)
                ind2 = tau2+kplus(2)+1;
                temp = A2(ind1:(ind1+k(1)-1), ind2:(ind2+k(2))-1,:);
                score(ind1,ind2) = norm(temp(:), 1);
            end
        end
        [temp,ind1] = max(score); [~,ind2] = max(temp);
        tau = [ind1(ind2) ind2]-kplus-1;
        A2 = circshift(A2,-tau);
        X2sol.X = circshift(X2sol.X,tau);
        X2sol.W = circshift(X2sol.W,tau);

        % Save phase 2 extras:
        if i == 1;  idx = 1;    else; idx = i;    end
        extras.phase2(idx).A = A2;
        extras.phase2(idx).X = X2sol.X;
        extras.phase2(idx).b = X2sol.b;
        extras.phase2(idx).info = info;
        if i == 1;  extras.phase2 = fliplr(extras.phase2);  end

        dispfun2(A2,X2sol.X);
        lambda = lambda*lam2fac;
        i = i+1;

    end
end

%% Finished: get the final A, X
if params.phase2
    Aout = A2(kplus(1)+(1:k(1)), kplus(2)+(1:k(2)), :);
    extras.normA = norm(Aout(:));
    Xout = circshift(X2sol.X,kplus) * norm(Aout(:));
    Aout = Aout/norm(Aout(:));
    bout = X2sol.b;
else
    Aout = A;
    extras.normA = norm(Aout(:));
    Xout = Xsol.X;
    bout = Xsol.b;
end

runtime = toc(starttime);
fprintf('\nDone! Runtime = %.2fs. \n\n', runtime);
end


function [ Aout, Xsol, extras ] = Asolve_Manopt( Y, Ain, lambda, Xsolve, varargin )
%ASolve_MANOPT     BD using Manopt solvers.
%   - Core usage:
%       [ Aout, Xsol, Stats ] = Asolve_Manopt( Y, Ain, lambda, Xsolve)
%
%   - Optional variables:
%       [ ... ] = Asolve_Manopt( ... , Xinit, Xpos, getbias, dispfun )
%

    load([fileparts(mfilename('fullpath')) '/../config/Asolve_config.mat']); %#ok<*LOAD>

    k = size(Ain);
    if (numel(k) > 2)
        n = k(3); k = k(1:2);
    else
        n = 1;
    end

    Ain = Ain/norm(Ain(:));

    %% Handle the extra variables:
    nvarargin = numel(varargin);
    if nvarargin > 4
        error('Too many input arguments.');
    end

    idx = 3;
    if nvarargin < idx || isempty(varargin{idx})
        getbias = false;
    else
        getbias = varargin{idx};
    end
    
    idx = 2;
    if nvarargin < idx || isempty(varargin{idx})
        xpos = false;
    else
        xpos = varargin{idx};
    end
    
    idx = 1;
    if nvarargin < idx || isempty(varargin{idx})
        if strcmp(Xsolve,'FISTA')
            xinit = Xsolve_FISTA(Y, Ain, lambda, mu, [], xpos);
        elseif strcmp(Xsolve,'pdNCG')
            xinit = Xsolve_pdNCG(Y, Ain, lambda, mu, [], xpos);
        end
    else
        xinit = varargin{idx};
    end

    idx = 4;
    if nvarargin < idx || isempty(varargin{idx})
        dispfun = @(a, X) 0;
    else
        dispfun = varargin{idx};
    end

    %% Set up the problem structure for Manopt and solve
    % The package containing supplement information for the cost, egrad,
    % and ehess functions:
    %{
    suppack.Y = Y;
    suppack.k = k;
    suppack.n = n;
    suppack.mu = mu;
    suppack.xinit = xinit;
    suppack.saveiterates = saveiterates;
    %}

    problem.M = spherefactory(prod(k)*n);
    problem.cost = @(a, store) costfun(a, store, Y, k, n, lambda, mu, xinit, xpos, getbias, Xsolve);
    problem.egrad = @(a, store) egradfun(a, store,  Y, k, n, lambda, mu, xinit, xpos, getbias, Xsolve);
    problem.ehess = @(a, u, store) ehessfun(a, u, store,  Y, k, n, lambda, mu, xinit, xpos, getbias, Xsolve);

    options.statsfun = @(problem, a, stats, store) statsfun( problem, a, stats, store, k, n, saveiterates, dispfun);
    %options.stopfun = @(problem, x, info, last) stopfun(problem, x, info, last, TRTOL);

    % Run Manopt solver:
    [Aout, extras.cost, info, extras.options] = ManoptSolver(problem, Ain(:), options);

    % Produce final output:
    Aout = reshape(Aout, [k n]);
    if saveiterates
        extras.Aiter = arrayfun(@(i) info(i).A, 1:numel(info), 'UniformOutput', false);
        niter = numel(extras.Aiter);
        extras.Aiter = cell2mat(reshape(extras.Aiter, [1 1 niter]));
        extras.Xiter = arrayfun(@(i) info(i).X, 1:numel(info), 'UniformOutput', false);
        extras.Xiter = cell2mat(reshape(extras.Xiter, [1 1 niter]));

        Xsol.X = extras.Xiter(:,:,end);
        Xsol.W = info(end).W;
        Xsol.b = info(end).b;
    else
        if strcmp(Xsolve,'FISTA')
            Xsol = Xsolve_FISTA( Y, Aout, lambda, mu, xinit, xpos, getbias);
        elseif strcmp(Xsolve,'pdNCG')
            Xsol = Xsolve_pdNCG(Y, Aout, lambda, mu, xinit, xpos, getbias);
        end
    end
end

function [ cost, store ] = costfun( a, store, Y, k, n, lambda, mu, xinit, xpos, getbias, Xsolve)
    if ~isfield(store, 'X')
        store = computeX( a, store, Y, k, n, lambda, mu, xinit, xpos, getbias, Xsolve);
    end

    cost = store.cost;
end

function [ egrad, store ] = egradfun( a, store, Y, k, n, lambda, mu, xinit, xpos, getbias, Xsolve)
    if ~isfield(store, 'X')
        store = computeX( a, store, Y, k, n, lambda, mu, xinit, xpos, getbias, Xsolve);
    end
    
    m = size(store.X);
    
    egrad = zeros(prod(k)*n,1);
    
    for i = 1:n
        idx = (i-1)*prod(k) + (1:prod(k));
        tmp = convfft2( store.X, convfft2( reshape(a(idx), k), store.X ) + store.b(i) - Y(:,:,i), 1, m+k-1, m);
        tmp = tmp(1:k(1), 1:k(2));
        egrad(idx) = tmp(:);
    end
end

function [ ehess, store ] = ehessfun( a, u, store, Y, k, n, lambda, mu, xinit, xpos, getbias, Xsolve)
    if ~isfield(store, 'X')
        store = computeX( a, store, Y, k, n, lambda, mu, xinit, xpos, getbias, Xsolve);
    end

    ehess = H_function( u, Y, reshape(a, [k n]), store.X, lambda, mu );
end

function [ store ] = computeX( a, store, Y, k, n, lambda, mu, xinit, xpos, getbias, Xsolve)
    % Updates the cache to store X*(A), and the active-set whenever a new
    % a new iteration by the trust-region method needs it.
    if strcmp(Xsolve,'FISTA')
        sol = Xsolve_FISTA( Y, reshape(a, [k n]), lambda, mu, xinit, xpos, getbias );
    elseif strcmp(Xsolve,'pdNCG')
        sol = Xsolve_pdNCG( Y, reshape(a, [k n]), lambda, mu, xinit, xpos, getbias);
    end
        
    store.X = sol.X;
    store.W = sol.W;
    store.b = sol.b;
    store.cost = sol.f;
end

function [ stats ] = statsfun( problem, a, stats, store, k, n, saveiterates, dispfun) %#ok<INUSL>
    if saveiterates
        stats.A = reshape(a, [k n]);
        stats.X = store.X;      % So X could be returned at the end.
        stats.W = store.W;
        stats.b = store.b;
    end
    dispfun(a, store.X);
end

function [ Xsol, info ] = Xsolve_FISTA( Y, A, lambda, mu, varargin )
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

    
    
    %% Iterate:    
    t=1; W = X; u = b;
    costs = NaN(MAXIT,2);
    doagain = true;  it = 0;  count = 0;
    while doagain
	it = it + 1;
        % Gradients and Hessians:
        grad_fW = zeros(m); grad_fu = zeros(n,1); R_A = zeros(m);
        for i = 1:n     % sum up
            Ri = convfft2(A(:,:,i), W) + u(i) - Y(:,:,i);
            grad_fW = grad_fW + convfft2( A(:,:,i), Ri, 1 );
            grad_fu(i) = sum(Ri(:));
            R_A = R_A + abs(fft2(A(:,:,i),m(1),m(2))).^2;
        end

        % FISTA update
        L = max(R_A(:));
        %size(W);
        X_ = g.prox(W - 1/L*grad_fW, lambda/L, xpos);
        %size(X_);
        t_ = (1+sqrt(1+4*t^2))/2;
        W = X_ + (t-1)/t_*(X_-X);
        if getbias
            b_ = u - grad_fu/(2*prod(m)*sqrt(n));  
            u = b_ + (t-1)/t_*(b_-b);
            b = b_;
        end
        X = X_; t = t_;

        %TODO Check conditions to repeat iteration:
        f = 0;
        for i = 1:n
            f = f + norm(convfft2(A(:,:,i), reshape(X, m)) + b(i) - Y(:,:,i), 'fro')^2/2;
        end
        costs(it,1) = f;
        costs(it,2) = g.cost(X, lambda);

        tmp = grad_fW;
        for i = 1:n
            %tmp(:,:,i) = tmp(:,:,i) + grad_fu(i); 
            tmp = tmp + grad_fu(i); %testing
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
    Xsol.f = sum(costs(it,:));
    info.numit = it;
    info.delta = delta;
    info.costs = costs(1:it,:);
end

function [ C ] = cconvfft2( A, B, varargin)
%CCONVFFT2	FFT implementation of 2D cyclic convolution
%   C = cconvfft(A, B)	convolves A and B using the larger size
%
%   C = cconvfft(A, B, N)   convolves A and B using size N
%
%   C = cconvfft(A, B, N, adj)  convolves A and B using size N, 
%   'adj == 'left'' leads cconvfft2 to convolve B with the adjoint kernel 
%   of A, and 'adj == 'right'' convolves A with the adjoint of B.
%
%   Both N and adj can be left empty.

    numvararg = numel(varargin);
    
    if numvararg > 2
        error('Too many input arguments.');
    end
    
    N = max(size(A), size(B));
    if numvararg >= 1 && ~isempty(varargin{1})
        N = varargin{1};
    end
    
    A_hat = fft2(A,N(1),N(2));
    B_hat = fft2(B,N(1),N(2));
    if numvararg >= 2 && ~isempty(varargin{2})
        switch varargin{2}
            case 'left'
                A_hat = conj(A_hat);
            case 'right'
                B_hat = conj(B_hat);
        end
    end
    
    C = ifft2( A_hat .* B_hat );

end

function [ Y ] = convfft2( A, B, varargin )
%CCONVFFT2	FFT implementation of 2D linear convolution
%   Only non-extended portion of convolution (size of B) kept.
%
%   Y = convfft2(A, B)	convolves A with B.
%   In 1D vector notation,      c = P*C{iota*a}*(rho*b), 
%   where rho extends B to size(A)+size(B)-1, and P shrinks the 
%   window observed from the output of Ca*rho*b to SIZE(B).
%
%   Y = convfft2(..., adj, tmpsz, outsz)  
%   applies the adjoint operation if ADJ==TRUE: rho'*C{iota*a}'*P'.
%   TMPSZ and OUTSZ sets the extension and output sizes manually.
%

    if numel(varargin) > 3
        error('Too many input arguments.');
    end
    if numel(varargin) >= 1 && ~isempty(varargin{1})
        adj = varargin{1};
    else
        adj = false;
    end
    if numel(varargin) >= 2 && ~isempty(varargin{2})
        tmpsz = varargin{2};
    else
        tmpsz = size(B) + size(A) - 1;
    end
    if numel(varargin) >= 3 && ~isempty(varargin{3})
        outsz = varargin{3};
    else
        outsz = size(B);
    end
    
    if adj == false
        Y = shrink(cconvfft2(A, extend(B,tmpsz)), outsz);
    else
        Y = extend(cconvfft2(A, shrink(B, tmpsz, 1), tmpsz, 'left'), outsz, 1);
    end
end



function [ out ] = extend( in, outsz, varargin )
%EXTEND     The extension map and its adjoint.
    if numel(varargin) > 1
        error('Too many input arguments.');
    end
    if numel(varargin) == 1 && ~isempty(varargin{1})
        adj = varargin{1};
    else
        adj = false;
    end
    
    if adj == false
        if sum(outsz < size(in))
            error('OUTSZ needs to be geq than SIZE(IN) if ''ADJ == FALSE''.');
        end
        out = zeros(outsz);
        out(1:size(in,1), 1:size(in,2)) = in;
    else
        if sum(outsz > size(in))
            error('OUTSZ needs to be leq than SIZE(IN) if ''ADJ == TRUE''.');
        end
        out = in(1:outsz(1), 1:outsz(2));
    end

end

function [ H_v ] = H_function( v, Y, A, X, lambda, mu )
%H_FUNCTION     Apply the Euclidean Hessian to a vector.
    load([fileparts(mfilename('fullpath')) '/../../config/Hfunction_config.mat']);
    
    m = size(Y); k = size(A); n = size(Y,3);
    m = m(1:2); k = k(1:2);
    tmpsz = m + k - 1;
    
    Haa_v = zeros([k n]);
    r = zeros([m n]);
    Hxa_v = zeros(m);
    for i = 1:n
        idx = (i-1)*prod(k) + (1:prod(k));
        vi = reshape(v(idx), k);
        Haa_v(:,:,i) = convfft2( X, convfft2(X, vi, 0, tmpsz, m), 1, tmpsz, k);
    
        r(:,:,i) = convfft2(A(:,:,i),X) - Y(:,:,i);
        Hxa_v = Hxa_v ...
            + convfft2( A(:,:,i), convfft2(X, vi, 0, tmpsz, m), 1) ...
            + Hxres( r(:,:,i), vi, k, m );
    end
    
    hesspendiag = lambda * mu^2*(mu^2 + X(:).^2).^(-3/2);    
    Hxxfun = @(u) Hxx_function(u, m, A, hesspendiag);
    pcgprecond = @(u) u./(1 + hesspendiag);
    [HxxInv_Hxa_v,~] = pcg(Hxxfun, Hxa_v(:), PCGTOL, PCGIT, pcgprecond);
    HxxInv_Hxa_v = reshape(HxxInv_Hxa_v, m);
    
    % Hax = Ikm'*(CA*CX' + (CA*CX-CY)*Pi);
    
    Hax_HxxInv_Hxa_v = zeros([k n]);
    for i = 1:n
        Hax_HxxInv_Hxa_v(:,:,i) = ...
            convfft2( X, convfft2(A(:,:,i), HxxInv_Hxa_v, 0, tmpsz, m), 1, tmpsz, k) ...
            + Hxres( r(:,:,i), HxxInv_Hxa_v, m, k );
    end
    
    H_v = Haa_v - Hax_HxxInv_Hxa_v;
    H_v = H_v(:);

end

function [ out ] = Hxres( res, in, insz, outsz )
    tmpsz = outsz + insz - 1;
    Vhat = fft2(extend(reshape(in, insz), tmpsz));
    reshat = fft2(shrink(res, tmpsz, 1));
    out = extend(ifft2(reshat.*conj(Vhat)), outsz, 1 );
end

function [ h ] = huber( mu )
    if nargin < 1 || isempty(mu)
        mu = 1e-6;
    end

    h.cost = @(x, lambda) cost(x, lambda, mu);
    h.prox = @(x, lambda, xpos) prox (x, lambda, mu, xpos);
    h.diffsubg = @(x, y, lambda, xpos) diffsubg(x, y, lambda, mu, xpos);
end

function [ hx ] = cost(x, lambda, mu)
    leq = abs(x) <= mu;
    hx = sum(x(leq).^2/2/mu) + sum(abs(x(~leq)) - mu/2);
    hx = lambda*hx;
end

function [ proxh ] = prox(x, lambda, mu, xpos)
    if nargin < 4 || isempty(xpos)
        xpos = false;
    end

    leq = abs(x) <= lambda + mu;
    proxh = zeros(size(x));
    proxh(leq) = x(leq)./(1+lambda/mu);
    proxh(~leq) = x(~leq) - lambda*sign(x(~leq));

    if xpos;  proxh = max(proxh,0);  end
end

function [ eps ] = diffsubg(x, y, lambda, mu, xpos)
%DIFFSUBG  Difference from y to the subgradient of the loss function at x.
%  In this case the Huber loss function is smooth so only need to subtract
%  from the gradient.

    leq = abs(x) <= mu;
    subg = zeros(size(x));

    subg(leq) = x(leq)/mu;
    subg(~leq) = sign(x(~leq));
    if xpos;  subg(x<0) = 0;  end

    eps = y - lambda*subg;
end


function [ out ] = Hxx_function( v, m, A, Hdiag )
    %% 1. Apply rho'*Ca'*P'*P*Ca*rho to v
    % Reshape and extend v
    k = [size(A,1) size(A,2)]; n = size(A,3);
    tmpsz = m + k - 1;
    rhov = extend(reshape(v, m), tmpsz);
    
    % For each slice, apply the rest of the operations
    tmp = zeros(m);
    delta = tmpsz - m;
    mask = false(tmpsz); 
    mask(floor(delta(1)/2)+(1:m(1)), floor(delta(2)/2)+(1:m(2))) = true;
    for i = 1:n
        % apply P'*P*Ca to rhov
        currslice = mask.*cconvfft2(A(:,:,i),rhov);
        
        % apply rho'*Ca' and save
        currslice = extend(cconvfft2(A(:,:,i),currslice,tmpsz,'left'),m,1);
        tmp = tmp + currslice;
    end
    out = tmp(:) + Hdiag.*v;
end

function [ Anorm ] = proj2oblique( A )
%PROJ2OBLIQUE   Normalize each slice to lie on the sphere.

    Anorm = NaN(size(A));
    for i = 1:size(A,3)
        Anorm(:,:,i) = A(:,:,i)/norm(A(:,:,i),'fro');
    end

end

function [ out ] = shrink( in, outsz, varargin )
%RESTRICT   Restriction of observation window and its adjoint.
    if numel(varargin) > 1
        error('Too many input arguments.');
    end
    if numel(varargin) == 1 && ~isempty(varargin{1})
        adj = varargin{1};
    else
        adj = false;
    end
    
    if adj == false
        if sum(outsz > size(in))
            error('OUTSZ needs to be leq than SIZE(IN) if ''ADJ == FALSE''.');
        end
        delta = size(in) - outsz;
        out = in(floor(delta(1)/2)+(1:outsz(1)), floor(delta(2)/2)+(1:outsz(2)));
    else
        if sum(outsz < size(in))
            error('OUTSZ needs to be geq than SIZE(IN) if ''ADJ == TRUE''.');
        end
        delta = outsz - size(in);
        out = zeros(outsz);
        out(floor(delta(1)/2)+(1:size(in,1)), floor(delta(2)/2)+(1:1:size(in,2))) = in;
    end
end

