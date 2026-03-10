function [ Aout, Xout, bout, extras ] = MTSBD_all_slice_modified( Y, k, params, dispfun, kernel_initialguess, Max_iteration, maxIT)
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
    %   Finally, two optional fields for the struct. These features are
    %   automatically disabled if the fields are not included or are empty:
    %
    %       Xsolve,   string     : Pick which Xsolve to use--'FISTA' or 'pdNCG.'
    %       xpos,     bool       : Constrain X to have nonnegative entries
    %       getbias,  bool       : Extract constant bias from observation.
    %                                                                                                                      
    %   k: n*2 matrix containing a list of n different kernel sizes [x,y]. 
    %   kernel_initialguess: cell array containing different intial guess of kernel 
    
    %% Start timing for the whole process
    total_starttime = tic;
    t_config_start = tic;
    
    %% Process input arguments
    if nargin < 4 || isempty(dispfun)
        dispfun = @(Y,A,X,k,kplus,idx) 0;
    end
    
    if isfield(params, 'lambda1_weighted') && ~isempty(params.lambda1_weighted)
        lambda1_weighted = params.lambda1_weighted;
    else
        lambda1_weighted = params.lambda1;
    end
    if isfield(params, 'lambda1_unweighted') && ~isempty(params.lambda1_unweighted)
        lambda1_unweighted = params.lambda1_unweighted;
    else
        lambda1_unweighted = params.lambda1;
    end
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
    
    if ~isfield(params, 'xinit') || isempty(params.xinit)
        xinit = [];
    else
        xinit = params.xinit;
    end

    if ~isfield(params, 'noise_var')
        error('params must contain noise_var for quality metrics');
    end
    noise_var = params.noise_var;

    slices = size(Y,3);
    spatial = size(Y,1);
    kernel_num = size(k,1);
    mu = 10^-6;

    if ~isfield(params, 'kernel_update_order') || isempty(params.kernel_update_order)
        kernel_update_order = 1:kernel_num;
    else
        kernel_update_order = params.kernel_update_order(:).';
        if numel(kernel_update_order) ~= kernel_num
            error('params.kernel_update_order must contain exactly one index per kernel.');
        end
        if any(kernel_update_order < 1) || any(kernel_update_order > kernel_num)
            error('params.kernel_update_order indices must be in [1, kernel_num].');
        end
        if numel(unique(kernel_update_order)) ~= kernel_num
            error('params.kernel_update_order must be a permutation without duplicates.');
        end
    end

    if ~isfield(params, 'slice_weights') || isempty(params.slice_weights)
        slice_weights = ones(slices, kernel_num);
    else
        sw = params.slice_weights;
        if isvector(sw)
            sw = sw(:);
            if numel(sw) ~= slices
                error('Vector params.slice_weights must have one value per slice.');
            end
            slice_weights = repmat(sw, [1, kernel_num]);
        else
            if ~isequal(size(sw), [slices, kernel_num])
                error('Matrix params.slice_weights must be [num_slices x num_kernels].');
            end
            slice_weights = sw;
        end

        % Validate each kernel has at least one trusted slice.
        for n = 1:kernel_num
            if any(slice_weights(:,n) < 0)
                error('slice_weights must be nonnegative.');
            end
            if sum(slice_weights(:,n)) <= 0
                error('Each kernel weight column must contain at least one positive entry.');
            end
        end
    end
    
    % Update the configuration file with the new max_iteration
    update_config('Xsolve_config.mat', 'MAXIT', Max_iteration, 'Xsolve_config_tunable.mat');
    update_config('Asolve_config.mat','options.maxiter', Max_iteration, 'Asolve_config_tunable.mat');
    t_config = toc(t_config_start);
    
    %% Phase I: Initialization and First Iteration with demixing applied
    fprintf('PHASE I: Initialization and First Iteration\n');
    A = kernel_initialguess;
    X_struct = struct();
    Xiter = zeros([spatial,spatial,kernel_num]);
    biter = zeros(slices, kernel_num);

    % use the initial guess if provided
    if ~isempty(xinit)
        for n = 1:kernel_num
            Xiter(:,:,n) = xinit{n}.X;
            biter(:,n) = xinit{n}.b;
        end
    end
    
    % Demixing blend schedule (depends on maxIT):
    % faint_factor(iter) transitions from 0.5 (first iteration) to 0.1 (last).
    faint_factor_first = 0.5;
    faint_factor_last = 0.1;
    
    kernel_time_total = 0;
    total_phase1_sum = 0;
    
    % Main iteration loop
    for iter = 1:maxIT
        iter_starttime = tic;
        
        % Compute Y_background with demixing approach

        Y_mainloop_perkernel = zeros(size(Y,1),size(Y,2),size(Y,3),kernel_num);

        for m = 1:kernel_num
            Y_mainloop_perkernel(:,:,:,m) = convfft3(A{m},Xiter(:,:,m)) + reshape(biter(:,m),1,1,size(biter,1));
        end

        % Update each kernel
        for ord_idx = 1:kernel_num
            t_miniloop_start = tic;
            n = kernel_update_order(ord_idx);
            Y_sum = sum(Y_mainloop_perkernel,4);
            % update Y_residual 
            Y_residual = Y - Y_sum;
            if maxIT > 1
                faint_factor = faint_factor_first + (faint_factor_last - faint_factor_first) * ((iter - 1) / (maxIT - 1));
            else
                faint_factor = faint_factor_first;
            end
            Yiter = Y_residual + (1-faint_factor)*Y_mainloop_perkernel(:,:,:,n) + faint_factor*Y_sum;
            dispfun1 = @(A, X) dispfun{n}(Y(:,:,1), A(:,:,1), X, k(n,:), []);
            
            % Initial X computation only on first iteration if no initial guess is provided
            if iter == 1
                if isempty(xinit)
                    X_struct.(['x',num2str(n)]) = Xsolve_FISTA_tunable(Y, A{n}, lambda1_weighted(n), mu, xinit, xpos, getbias, slice_weights(:,n));
                else
                    X_struct.(['x',num2str(n)]) = Xsolve_FISTA_tunable(Y, A{n}, lambda1_weighted(n), mu, xinit{n}, xpos, getbias, slice_weights(:,n));
                end
            end
            
            [A{n}, X_struct.(['x',num2str(n)]), info] = Asolve_Manopt_tunable(Yiter, A{n}, lambda1_unweighted(n), Xsolve, X_struct.(['x',num2str(n)]), xpos, getbias, dispfun1, slice_weights(:,n));
            
            Xiter(:,:,n) = X_struct.(['x',num2str(n)]).X;
            biter(:,n) = X_struct.(['x',num2str(n)]).b;
            % update Y_mainloop_perkernel for this kernel type 
            Y_mainloop_perkernel(:,:,:,n)=convfft3(A{n},Xiter(:,:,n))+reshape(biter(:,n),1,1,size(biter,1));
            
            ktime = toc(t_miniloop_start);
            kernel_time_total = kernel_time_total + ktime;
            fprintf('kernel %d finished, runtime = %.2fs\n', n, ktime);
        end
        % Recompute residual after all kernel updates so metrics use
        % the final state of this outer iteration (not a stale inner-loop value).
        Y_sum = sum(Y_mainloop_perkernel,4);
        Y_residual = Y - Y_sum;
        residual_var = reshape(var(Y_residual, 0, [1,2]),[1,size(Y,3)]);
        % Calculate quality metric
        quality_metric = mean(noise_var ./ residual_var);
        extras.phase1.quality_metrics(iter) = quality_metric;
        extras.phase1.residuals(:,:,:,iter) = Y_residual;
        iter_runtime = toc(iter_starttime);
        total_phase1_sum = total_phase1_sum + iter_runtime;
        fprintf('Iteration %d: Runtime = %.2fs\n', iter, iter_runtime);
    end
    
    % Store results
    extras.phase1.Aout = A;
    extras.phase1.Xout = Xiter;
    extras.phase1.biter = biter;
    
    %% PHASE II: Lift the sphere and do lambda continuation
    if params.phase2
        fprintf('\n\nPHASE II: \n=========\n');
        k3 = k + 2*kplus;
    
        A2 = cell(1, kernel_num);
        X2_struct = struct();
        for n = 1:kernel_num
            A2{n} = zeros(k3(n,:));
            A2{n}(kplus(n,1)+(1:k(n,1)), kplus(n,2)+(1:k(n,2))) = A{n};
            X2_struct.(['x',num2str(n)]) = X_struct.(['x',num2str(n)]);
        end
    
        lambda = lambda1_unweighted;
        lam2fac = (lambda2./lambda1_unweighted).^(1/nrefine);
        
        % Initialize Phase II metrics
        extras.phase2.activation_metrics = zeros(nrefine + 1, kernel_num);
        extras.phase2.kernel_quality_factors = zeros(nrefine + 1, kernel_num);
        
        for i = 1:nrefine + 1
            fprintf('lambda iteration %d/%d: \n', i, nrefine + 1);
            
            % Build per-kernel reconstruction cache for this refinement.
            Y_phase2_perkernel = zeros([size(Y), kernel_num]);
            for m = 1:kernel_num
                Y_phase2_perkernel(:,:,:,m) = convfft2(A2{m}, X2_struct.(['x',num2str(m)]).X) ...
                    + reshape(X2_struct.(['x',num2str(m)]).b, 1, 1, []);
            end
    
            for ord_idx = 1:kernel_num
                n = kernel_update_order(ord_idx);
                fprintf('Processing kernel %d, lambda = %.1e: \n', n, lambda(n));
                % Calculate Yiter without demixing
                Y_sum = sum(Y_phase2_perkernel, 4);
                Y_residual = Y - Y_sum;
                Yiter = Y_residual + Y_phase2_perkernel(:,:,:,n);
                
                dispfun2 = @(A, X) dispfun{n}(Y, A, X, k3(n,:), kplus(n,:));
                [A2{n}, X2_struct.(['x',num2str(n)]), info] = Asolve_Manopt_tunable(Yiter, A2{n}, lambda(n), Xsolve, X2_struct.(['x',num2str(n)]), xpos, getbias, dispfun2, slice_weights(:,n));
                
                % Attempt to 'unshift" the a and x
                score = zeros(2*kplus(n,1)+1, 2*kplus(n,2)+1);
                for tau1 = -kplus(n,1):kplus(n,1)
                    ind1 = tau1+kplus(n,1)+1;
                    for tau2 = -kplus(n,2):kplus(n,2)
                        ind2 = tau2+kplus(n,2)+1;
                        temp = A2{n}(ind1:(ind1+k(n,1)-1), ind2:(ind2+k(n,2)-1));
                        score(ind1,ind2) = norm(temp(:), 1);
                    end
                end
                [temp,ind1] = max(score); [~,ind2] = max(temp);
                tau = [ind1(ind2) ind2]-kplus(n,:)-1;
                A2{n} = circshift(A2{n},-tau);
                X2_struct.(['x',num2str(n)]).X = circshift(X2_struct.(['x',num2str(n)]).X,tau);
                X2_struct.(['x',num2str(n)]).W = circshift(X2_struct.(['x',num2str(n)]).W,tau);

                % Refresh only the updated kernel contribution.
                Y_phase2_perkernel(:,:,:,n) = convfft2(A2{n}, X2_struct.(['x',num2str(n)]).X) ...
                    + reshape(X2_struct.(['x',num2str(n)]).b, 1, 1, []);
    
                % Save phase 2 extras
                extras.phase2.A{n} = A2{n};
                extras.phase2.X{n} = X2_struct.(['x',num2str(n)]).X;
                extras.phase2.b(n) = X2_struct.(['x',num2str(n)]).b;
                extras.phase2.info{n} = info;
            end
            
            % Evaluate metrics for this refinement
            X2_combined = zeros(size(Y,1), size(Y,2), kernel_num);
            A2_central = cell(1, kernel_num);
            for n = 1:kernel_num
                X2_combined(:,:,n) = X2_struct.(['x',num2str(n)]).X;
                A2_central{n} = A2{n}(kplus(n,1)+(1:k(n,1)), kplus(n,2)+(1:k(n,2)));
            end
    
            [activation_similarity, kernel_similarity] = computeQualityMetrics(X0, X2_combined, A0, A2_central, k3);
            extras.phase2.activation_metrics(i,:) = activation_similarity;
            extras.phase2.kernel_quality_factors(i,:) = kernel_similarity;
            
            % Print results
            fprintf('Refinement %d Metrics:\n', i);
            fprintf('Activation Quality:\n');
            for n = 1:kernel_num
                fprintf('Kernel %d - Similarity: %.3f\n', n, activation_similarity(n));
            end
            fprintf('Kernel Quality Factors:\n');
            for n = 1:kernel_num
                fprintf('Kernel %d: %.3f\n', n, kernel_similarity(n));
            end
            
            lambda = lambda .* lam2fac;
        end
    end
    
    %% Finished: get the final A, X
    if params.phase2
        Aout = cell(1, kernel_num);
        Xout = zeros(size(Y,1), size(Y,2), kernel_num);
        bout = zeros(kernel_num, 1);
        extras.normA = zeros(kernel_num, 1);
        
        for n = 1:kernel_num
            Aout{n} = A2{n}(kplus(n,1)+(1:k(n,1)), kplus(n,2)+(1:k(n,2)));
            extras.normA(n) = norm(Aout{n}(:));
            Xout(:,:,n) = circshift(X2_struct.(['x',num2str(n)]).X, kplus(n,:)) * extras.normA(n);
            Aout{n} = Aout{n} / extras.normA(n);
            bout(n) = X2_struct.(['x',num2str(n)]).b;
        end
    else
        Aout = A;
        extras.normA = zeros(kernel_num, 1);
        for n = 1:kernel_num
            extras.normA(n) = norm(Aout{n}(:));
        end
        Xout = Xiter;
        bout = biter;
    end
    
    %% Final timing
    total_runtime = toc(total_starttime);
    fprintf('\nTotal Runtime = %.2fs\n\n', total_runtime);
    
    % Time breakdown (Phase I only; Phase II if run is included in total_runtime)
    phase1_other = total_phase1_sum - kernel_time_total;
    phase2_and_rest = total_runtime - t_config - total_phase1_sum;
    fprintf('Time breakdown (%% of total):\n');
    fprintf('  Config update:           %.2fs (%5.1f%%)\n', t_config, 100*t_config/max(total_runtime,1e-6));
    fprintf('  Phase I kernel updates:  %.2fs (%5.1f%%)  [Asolve_Manopt + inner Xsolve]\n', kernel_time_total, 100*kernel_time_total/max(total_runtime,1e-6));
    fprintf('  Phase I other:           %.2fs (%5.1f%%)  [demixing, Y_residual, quality metric]\n', phase1_other, 100*phase1_other/max(total_runtime,1e-6));
    if phase2_and_rest > 0.01
        fprintf('  Phase II + rest:         %.2fs (%5.1f%%)\n', phase2_and_rest, 100*phase2_and_rest/max(total_runtime,1e-6));
    end
    fprintf('  Total:                   %.2fs (100%%)\n', total_runtime);
    
    % Store runtime in extras
    extras.runtime = total_runtime;
end