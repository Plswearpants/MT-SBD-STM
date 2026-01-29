function [ Aout, Xsol, extras ] = Asolve_Manopt_ALL( Y, Ain, lambda, Xsolve, varargin )
%ASolve_MANOPT     BD using Manopt solvers.
%   - Core usage:
%       [ Aout, Xsol, Stats ] = Asolve_Manopt( Y, Ain, lambda, Xsolve)
%
%   - Optional variables:
%       [ ... ] = Asolve_Manopt( ... , Xinit, Xpos, getbias, dispfun )
%

    load([fileparts(mfilename('fullpath')) '/../../config/Asolve_config.mat']); %#ok<*LOAD>

    k = size(Ain);
    if (numel(k) > 3)
        t = k(1); % number of types
        n = k(4); % number of slices
        k = k(2:3); % kernel size
    elseif (numel(k) > 2)
        t = 1; % single type
        n = k(3); % number of slices
        k = k(1:2); % kernel size
        Ain = reshape(Ain, [t, k, n]); % reshape to match new dimension format
    else
        t = 1; % single type
        n = 1; % single slice
        k = k(1:2); % kernel size
        Ain = reshape(Ain, [t, k, n]); % reshape to match new dimension format
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

    problem.M = spherefactory(t*prod(k)*n);
    problem.cost = @(a, store) costfun(a, store, Y, t, k, n, lambda, mu, xinit, xpos, getbias, Xsolve);
    problem.egrad = @(a, store) egradfun(a, store, Y, t, k, n, lambda, mu, xinit, xpos, getbias, Xsolve);
    problem.ehess = @(a, u, store) ehessfun(a, u, store, Y, t, k, n, lambda, mu, xinit, xpos, getbias, Xsolve);

    options.statsfun = @(problem, a, stats, store) statsfun(problem, a, stats, store, t, k, n, saveiterates, dispfun);
    %options.stopfun = @(problem, x, info, last) stopfun(problem, x, info, last, TRTOL);

    % Run Manopt solver:
    [Aout, extras.cost, info, extras.options] = ManoptSolver(problem, Ain(:), options);

    % Produce final output:
    Aout = reshape(Aout, [t, k, n]);
    if saveiterates
        extras.Aiter = arrayfun(@(i) info(i).A, 1:numel(info), 'UniformOutput', false);
        niter = numel(extras.Aiter);
        extras.Aiter = cell2mat(reshape(extras.Aiter, [1 1 1 niter]));
        extras.Xiter = arrayfun(@(i) info(i).X, 1:numel(info), 'UniformOutput', false);
        extras.Xiter = cell2mat(reshape(extras.Xiter, [1 1 1 niter]));

        Xsol.X = extras.Xiter(:,:,:,end);
        Xsol.W = info(end).W;
        Xsol.b = info(end).b;
    else
        if strcmp(Xsolve,'FISTA')
            Xsol = Xsolve_FISTA(Y, Aout, lambda, mu, xinit, xpos, getbias);
        elseif strcmp(Xsolve,'pdNCG')
            Xsol = Xsolve_pdNCG(Y, Aout, lambda, mu, xinit, xpos, getbias);
        end
    end
end

function [ cost, store ] = costfun(a, store, Y, t, k, n, lambda, mu, xinit, xpos, getbias, Xsolve)
    if ~isfield(store, 'X')
        store = computeX(a, store, Y, t, k, n, lambda, mu, xinit, xpos, getbias, Xsolve);
    end

    cost = store.cost;
end

function [ egrad, store ] = egradfun(a, store, Y, t, k, n, lambda, mu, xinit, xpos, getbias, Xsolve)
    if ~isfield(store, 'X')
        store = computeX(a, store, Y, t, k, n, lambda, mu, xinit, xpos, getbias, Xsolve);
    end
    
    m = size(store.X);
    m = m(2:3); % spatial dimensions
    
    egrad = zeros(t*prod(k)*n, 1);
    
    for i = 1:n
        for j = 1:t
            idx = ((i-1)*t + (j-1))*prod(k) + (1:prod(k));
            
            % Compute residual: Y - sum_over_types(conv(A, X))
            residual = Y(:,:,i);
            for typ = 1:t
                typ_a_idx = ((i-1)*t + (typ-1))*prod(k) + (1:prod(k));
                residual = residual - convfft2(reshape(a(typ_a_idx), k), store.X(typ,:,:,i));
            end
            
            if getbias
                residual = residual - store.b(i);
            end
            
            % Compute gradient for this type and slice
            tmp = convfft2(store.X(j,:,:,i), residual, 1, m+k-1, m);
            tmp = tmp(1:k(1), 1:k(2));
            egrad(idx) = -tmp(:); % Negative because residual is Y - conv(A,X)
        end
    end
end

function [ ehess, store ] = ehessfun(a, u, store, Y, t, k, n, lambda, mu, xinit, xpos, getbias, Xsolve)
    if ~isfield(store, 'X')
        store = computeX(a, store, Y, t, k, n, lambda, mu, xinit, xpos, getbias, Xsolve);
    end

    ehess = H_function(u, Y, reshape(a, [t, k, n]), store.X, lambda, mu);
end

function [ store ] = computeX(a, store, Y, t, k, n, lambda, mu, xinit, xpos, getbias, Xsolve)
    % Updates the cache to store X*(A), and the active-set whenever a new
    % iteration by the trust-region method needs it.
    if strcmp(Xsolve,'FISTA')
        sol = Xsolve_FISTA(Y, reshape(a, [t, k, n]), lambda, mu, xinit, xpos, getbias);
    elseif strcmp(Xsolve,'pdNCG')
        sol = Xsolve_pdNCG(Y, reshape(a, [t, k, n]), lambda, mu, xinit, xpos, getbias);
    end
        
    store.X = sol.X;
    store.W = sol.W;
    store.b = sol.b;
    store.cost = sol.f;
end

function [ stats ] = statsfun(problem, a, stats, store, t, k, n, saveiterates, dispfun) %#ok<INUSL>
    if saveiterates
        stats.A = reshape(a, [t, k, n]);
        stats.X = store.X;      % So X could be returned at the end.
        stats.W = store.W;
        stats.b = store.b;
    end
    dispfun(a, store.X);
end