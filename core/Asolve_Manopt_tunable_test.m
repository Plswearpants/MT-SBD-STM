function [ Aout, Xsol, extras ] = Asolve_Manopt_tunable_test( Y, Ain, lambda, Xsolve, X_all_flat, gamma, varargin)
    %ASolve_MANOPT_TEST     BD using Manopt solvers with sequential cross-correlation regularizer.
    %   - Core usage:
    %       [ Aout, Xsol, Stats ] = Asolve_Manopt_tunable_test( Y, Ain, lambda, Xsolve, X_all_flat, gamma)
    %
    %   - Optional variables:
    %       [ ... ] = Asolve_Manopt_tunable_test( ... , Xinit, Xpos, getbias, dispfun )
    %

    load([fileparts(mfilename('fullpath')) '/../examples/Asolve_config_tunable.mat']); %#ok<*LOAD>
    
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

    % Validate X_all_flat
    if isempty(X_all_flat) || ~isnumeric(X_all_flat) || ~ismatrix(X_all_flat)
        error('X_all_flat must be a non-empty 2D matrix');
    end

    % Validate gamma
    if ~isnumeric(gamma) || ~isscalar(gamma) || gamma < 0
        error('gamma must be a non-negative scalar');
    end

    idx = 2;
    if nvarargin < idx || isempty(varargin{idx})
        xpos = false;
    else
        xpos = varargin{idx};
    end
    
    idx = 3;
    if nvarargin < idx || isempty(varargin{idx})
        getbias = false;
    else
        getbias = varargin{idx};
    end
    
    idx = 1;
    if nvarargin < idx || isempty(varargin{idx})
        if strcmp(Xsolve,'FISTA_test')
            xinit = Xsolve_FISTA_test(Y, Ain, lambda, mu, X_all_flat, gamma, [], xpos, getbias);
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
    problem.M = spherefactory(prod(k)*n);
    problem.cost = @(a, store) costfun(a, store, Y, k, n, lambda, mu, xinit, xpos, getbias, Xsolve, gamma, X_all_flat);
    problem.egrad = @(a, store) egradfun(a, store, Y, k, n, lambda, mu, xinit, xpos, getbias, Xsolve, gamma, X_all_flat);
    problem.ehess = @(a, u, store) ehessfun(a, u, store, Y, k, n, lambda, mu, xinit, xpos, getbias, Xsolve, gamma, X_all_flat);

    options.statsfun = @(problem, a, stats, store) statsfun(problem, a, stats, store, k, n, saveiterates, dispfun);
    %options.stopfun = @(problem, x, info, last) stopfun(problem, x, info, last, TRTOL);

    % Run Manopt solver:
    [Aout, extras.cost, info, extras.options] = ManoptSolver(problem, Ain(:), options);
    extras.info = info;
    
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
        if strcmp(Xsolve,'FISTA_test')
            Xsol = Xsolve_FISTA_test(Y, Aout, lambda, mu, X_all_flat, gamma, xinit, xpos, getbias);
        elseif strcmp(Xsolve,'pdNCG')
            Xsol = Xsolve_pdNCG(Y, Aout, lambda, mu, xinit, xpos, getbias);
        end
    end
end

function [ cost, store ] = costfun(a, store, Y, k, n, lambda, mu, xinit, xpos, getbias, Xsolve, gamma, X_all_flat)
    if ~isfield(store, 'X')
        store = computeX(a, store, Y, k, n, lambda, mu, xinit, xpos, getbias, Xsolve, gamma, X_all_flat);
    end

    cost = store.cost;
end

function [ egrad, store ] = egradfun(a, store, Y, k, n, lambda, mu, xinit, xpos, getbias, Xsolve, gamma, X_all_flat)
    if ~isfield(store, 'X')
        store = computeX(a, store, Y, k, n, lambda, mu, xinit, xpos, getbias, Xsolve, gamma, X_all_flat);
    end
    
    m = size(store.X);
    
    egrad = zeros(prod(k)*n,1);
    
    for i = 1:n
        idx = (i-1)*prod(k) + (1:prod(k));
        tmp = convfft2(store.X, convfft2(reshape(a(idx), k), store.X) + store.b(i) - Y(:,:,i), 1, m+k-1, m);
        tmp = tmp(1:k(1), 1:k(2));
        egrad(idx) = tmp(:);
    end
end

function [ ehess, store ] = ehessfun(a, u, store, Y, k, n, lambda, mu, xinit, xpos, getbias, Xsolve, gamma, X_all_flat)
    if ~isfield(store, 'X')
        store = computeX(a, store, Y, k, n, lambda, mu, xinit, xpos, getbias, Xsolve, gamma, X_all_flat);
    end

    ehess = H_function(u, Y, reshape(a, [k n]), store.X, lambda, mu);
end

function [ store ] = computeX(a, store, Y, k, n, lambda, mu, xinit, xpos, getbias, Xsolve, gamma, X_all_flat)
    % Updates the cache to store X*(A), and the active-set whenever a new
    % iteration by the trust-region method needs it.
    if strcmp(Xsolve,'FISTA_test')
        [Xsol, info] = Xsolve_FISTA_test(Y, reshape(a, [k n]), lambda, mu, X_all_flat, gamma, xinit ,xpos, getbias);
        store.X = Xsol.X;
        store.W = Xsol.W;
        store.b = Xsol.b;
        store.cost = sum(info.costs(end,:));  % Get the total cost including cross-correlation
    elseif strcmp(Xsolve,'pdNCG')
        sol = Xsolve_pdNCG(Y, reshape(a, [k n]), lambda, mu, xinit, xpos, getbias);
        store.X = sol.X;
        store.W = sol.W;
        store.b = sol.b;
        store.cost = sol.f;
    end
end

function [ stats ] = statsfun(problem, a, stats, store, k, n, saveiterates, dispfun) %#ok<INUSL>
    if saveiterates
        stats.A = reshape(a, [k n]);
        stats.X = store.X;      % So X could be returned at the end.
        stats.W = store.W;
        stats.b = store.b;
    end
    
    dispfun(a, store.X);
end 