function [ Aout, Xsol, extras ] = Asolve_Manopt_tunable_test( Y, Ain, lambda, Xsolve, varargin)
    %ASolve_MANOPT_TEST     BD using Manopt solvers with sequential cross-correlation regularizer.
    %   - Core usage:
    %       [ Aout, Xsol, Stats ] = Asolve_Manopt_tunable_test( Y, Ain, lambda, Xsolve, max_iteration)
    %
    %   - Optional variables:
    %       [ ... ] = Asolve_Manopt_tunable_test( ... , Xinit, Xpos, getbias, gamma, current_idx, dispfun )
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
    if nvarargin > 6
        error('Too many input arguments.');
    end

    idx = 1;
    if nvarargin < idx || isempty(varargin{idx})
        if strcmp(Xsolve,'FISTA_test')
            xinit = Xsolve_FISTA_test(Y, Ain, lambda, mu, [], xpos, getbias, gamma, current_idx);
        elseif strcmp(Xsolve,'pdNCG')
            xinit = Xsolve_pdNCG(Y, Ain, lambda, mu, [], xpos);
        end
    else
        xinit = varargin{idx};
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

    idx = 4;
    if nvarargin < idx || isempty(varargin{idx})
        gamma = 1e-3;  % Default cross-correlation regularization parameter
    else
        gamma = varargin{idx};
    end

    idx = 5;
    if nvarargin < idx || isempty(varargin{idx})
        current_idx = 1;  % Default to first activation
    else
        current_idx = varargin{idx};
    end
    
    idx = 6;
    if nvarargin < idx || isempty(varargin{idx})
        dispfun = @(a, X) 0;
    else
        dispfun = varargin{idx};
    end

    %% Set up the problem structure for Manopt and solve
    problem.M = spherefactory(prod(k)*n);
    problem.cost = @(a, store) costfun(a, store, Y, k, n, lambda, mu, xinit, xpos, getbias, Xsolve, gamma, current_idx);
    problem.egrad = @(a, store) egradfun(a, store, Y, k, n, lambda, mu, xinit, xpos, getbias, Xsolve, gamma, current_idx);
    problem.ehess = @(a, u, store) ehessfun(a, u, store, Y, k, n, lambda, mu, xinit, xpos, getbias, Xsolve, gamma, current_idx);

    options.statsfun = @(problem, a, stats, store) statsfun(problem, a, stats, store, k, n, saveiterates, dispfun);
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
        if strcmp(Xsolve,'FISTA_test')
            Xsol = Xsolve_FISTA_test(Y, Aout, lambda, mu, xinit, xpos, getbias, gamma, current_idx);
        elseif strcmp(Xsolve,'pdNCG')
            Xsol = Xsolve_pdNCG(Y, Aout, lambda, mu, xinit, xpos, getbias);
        end
    end
end

function [ cost, store ] = costfun(a, store, Y, k, n, lambda, mu, xinit, xpos, getbias, Xsolve, gamma, current_idx)
    if ~isfield(store, 'X')
        store = computeX(a, store, Y, k, n, lambda, mu, xinit, xpos, getbias, Xsolve, gamma, current_idx);
    end

    cost = store.cost;
end

function [ egrad, store ] = egradfun(a, store, Y, k, n, lambda, mu, xinit, xpos, getbias, Xsolve, gamma, current_idx)
    if ~isfield(store, 'X')
        store = computeX(a, store, Y, k, n, lambda, mu, xinit, xpos, getbias, Xsolve, gamma, current_idx);
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

function [ ehess, store ] = ehessfun(a, u, store, Y, k, n, lambda, mu, xinit, xpos, getbias, Xsolve, gamma, current_idx)
    if ~isfield(store, 'X')
        store = computeX(a, store, Y, k, n, lambda, mu, xinit, xpos, getbias, Xsolve, gamma, current_idx);
    end

    ehess = H_function(u, Y, reshape(a, [k n]), store.X, lambda, mu);
end

function [ store ] = computeX(a, store, Y, k, n, lambda, mu, xinit, xpos, getbias, Xsolve, gamma, current_idx)
    % Updates the cache to store X*(A), and the active-set whenever a new
    % iteration by the trust-region method needs it.
    if strcmp(Xsolve,'FISTA_test')
        [Xsol, info] = Xsolve_FISTA_test(Y, reshape(a, [k n]), lambda, mu, xinit, xpos, getbias, gamma, current_idx);
        store.X = Xsol.X;
        store.W = Xsol.W;
        store.b = Xsol.b;
        store.cost = info.costs(end,1);  % Get the final cost from info
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