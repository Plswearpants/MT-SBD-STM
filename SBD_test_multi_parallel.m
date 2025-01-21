function [ Aout, Xout, bout, extras ] = SBD_test_multi_parallel( Y, k, fixed_params, kernel_initialguess, dataset_idx, parameter_combinations, param_idx, maxIT)
    %SBD Summary of this function goes here
    %
    % Inputs:
    %   Y: Input data
    %   k: Kernel sizes [n x 2]
    %   fixed_params: Structure containing fixed parameters
    %   kernel_initialguess: Initial kernel guess
    %   dataset_idx: Index of the current dataset
    %   parameter_combinations: Matrix of parameter combinations
    %   param_idx: Index of current parameter combination
    %   maxIT: Maximum number of iterations
    %
    %   PARAMS STRUCT:
    %   ===============
    %   The options struct should include the fields:
    %       lambda1,  float > 0  : regularization parameter for Phase I (tuned)
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
    
    %% Input validation
    if isempty(dataset_idx)
        error('dataset_idx must not be empty');
    end
    if ~isnumeric(dataset_idx) || ~isscalar(dataset_idx)
        error('dataset_idx must be a scalar number');
    end
    
    %% Start timing for the whole process
    total_starttime = tic;
    
    % Get project root directory and config paths
    root_dir = fileparts(mfilename('fullpath'));
    config_dir = fullfile(root_dir, 'config', 'worker_configs', sprintf('config_param_%d', param_idx));
    
    try
        %% Process input arguments
        % Extract parameters from parameter_combinations
        lambda1 = zeros(1,size(k,1));
        for n = 1:size(k,1)
            lambda1(n)=parameter_combinations(param_idx,1);
        end
 
        mini_loop = parameter_combinations(param_idx,2); % Second column is mini_loop
        kernel_num = size(k,1);
        mu = 10^-6;
        
        % Optional parameters
        if isfield(fixed_params, 'Xsolve')
            Xsolve = fixed_params.Xsolve;
        else
            Xsolve = 'FISTA';
        end
        
        if isfield(fixed_params, 'xpos')
            xpos = fixed_params.xpos;
        else
            xpos = false;
        end
        
        if isfield(fixed_params, 'getbias')
            getbias = fixed_params.getbias;
        else
            getbias = false;
        end
        
        % Extract X0 and A0 for quality metrics
        if ~isfield(fixed_params, 'X0') || ~isfield(fixed_params, 'A0')
            error('params must contain X0 and A0 for quality metrics');
        end
        X0 = fixed_params.X0;
        A0 = fixed_params.A0;
        
        %% Phase I: Initialization and First Iteration
        
        fprintf('Param_idx=%d, Using lambda1=%.3e, mini_loop=%d\n', ...
            param_idx, lambda1(1), round(mini_loop));
        fprintf('PHASE I: Initialization and First Iteration\n');
        

        % Add demixing factor
        faint_factor = 1;

        A = kernel_initialguess;
        X_struct = struct();
        Xiter = zeros([size(Y),kernel_num]);
        biter = zeros(kernel_num,1);
        
        % Initialize metrics
        extras.phase1.activation_metrics = zeros(maxIT, kernel_num);
        extras.phase1.kernel_quality_factors = zeros(maxIT, kernel_num);
        
        % Store parameter values in extras
        extras.parameters.lambda1 = lambda1;
        extras.parameters.mini_loop = mini_loop;
        extras.parameters.param_idx = param_idx;
        extras.parameters.dataset_idx = dataset_idx;
        
        % Main iteration loop
        for iter = 1:maxIT
            iter_starttime = tic;
            
            % Compute Y_background (changed to demixing approach)
            Y_sum = zeros(size(Y));
            for m = 1:kernel_num
                if iter > 1
                    Y_sum = Y_sum + convfft2(A{m}, Xiter(:,:,m));
                end
            end
            Y_residual = Y - Y_sum;
            
            % Update each kernel
            for n = 1:kernel_num
                % Calculate Yiter for this kernel (changed to demixing approach)
                Yiter = Y_residual + (iter > 1) * (1-1/(faint_factor*iter+1))*convfft2(A{n}, Xiter(:,:,n)) + (1/(faint_factor*iter+1))*Y_sum;
                
                if iter == 1
                    % Initial X computation with parallel version
                    X_struct.(['x',num2str(n)]) = Xsolve_FISTA_parallel(Y, A{n}, ...
                        lambda1(n), mu, [], xpos, getbias, config_dir);
                end
                
                % Use parallel version of Asolve
                [A{n}, X_struct.(['x',num2str(n)]), info] = Asolve_Manopt_parallel(Yiter, ...
                    A{n}, lambda1(n), Xsolve, X_struct.(['x',num2str(n)]), xpos, ...
                    getbias, config_dir);
                
                Xiter(:,:,n) = X_struct.(['x',num2str(n)]).X;
                biter(n) = X_struct.(['x',num2str(n)]).b;
            end
            
            % Add input validation
            if isempty(X0) || isempty(Xiter) || isempty(A0) || isempty(A) || isempty(k)
                error('One or more inputs to computeQualityMetrics are empty');
            end
            
            % Pass dataset index and parameter index to quality metrics
            [activation_similarity, kernel_similarity] = computeQualityMetrics(...
                X0, Xiter, A0, A, k, ...
                extras.parameters.dataset_idx, extras.parameters.param_idx);
            extras.phase1.activation_metrics(iter,:) = activation_similarity;
            extras.phase1.kernel_quality_factors(iter,:) = kernel_similarity;
            
            iter_runtime = toc(iter_starttime);
            fprintf('Iteration %d: Runtime = %.2fs\n', iter, iter_runtime);
        end
        
        % Store results
        extras.phase1.Aout = A;
        extras.phase1.Xout = Xiter;
        extras.phase1.biter = biter;
        
        %% PHASE II: Lift the sphere and do lambda continuation
        if fixed_params.phase2
            fprintf('\n\nPHASE II: \n=========\n');
            k3 = k + 2*kplus;
        
            A2 = cell(1, kernel_num);
            X2_struct = struct();
            for n = 1:kernel_num
                A2{n} = zeros(k3(n,:));
                A2{n}(kplus(n,1)+(1:k(n,1)), kplus(n,2)+(1:k(n,2))) = A{n};
                X2_struct.(['x',num2str(n)]) = X_struct.(['x',num2str(n)]);
            end
        
            lambda = lambda1;
            lam2fac = (lambda2./lambda1).^(1/nrefine);
            
            % Initialize Phase II metrics
            extras.phase2.activation_metrics = zeros(nrefine + 1, kernel_num);
            extras.phase2.kernel_quality_factors = zeros(nrefine + 1, kernel_num);
            
            for i = 1:nrefine + 1
                fprintf('lambda iteration %d/%d: \n', i, nrefine + 1);
                
                % Standard Y_residual calculation (no demixing)
                for m = 1:kernel_num
                    Y_sum = Y_sum + convfft2(A2{m}, X2_struct.(['x',num2str(m)]).X);
                end
                Y_residual = Y - Y_sum;
        
                for n = 1:kernel_num
                    fprintf('Processing kernel %d, lambda = %.1e: \n', n, lambda(n));
                    % Calculate Yiter without demixing
                    Yiter = Y_residual + convfft2(A2{n}, X2_struct.(['x',num2str(n)]).X);
                    
                    % Use parallel version of Asolve
                    [A2{n}, X2_struct.(['x',num2str(n)]), info] = Asolve_Manopt_parallel(Yiter, ...
                        A2{n}, lambda(n), Xsolve, X2_struct.(['x',num2str(n)]), xpos, ...
                        getbias, config_dir);
                    
                    % Attempt to "unshift" the a and x
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
        if fixed_params.phase2
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
        
        % Store runtime in extras
        extras.runtime = total_runtime;
        
        % When computing quality metrics, ensure X0 is properly formatted
        if isfield(fixed_params, 'X0')
            X0_eval = fixed_params.X0;  % This should be [height × width × num_kernels]
            A0_eval = fixed_params.A0;  % This should be cell array of kernels
        else
            X0_eval = [];
            A0_eval = [];
        end

        % Compute quality metrics with proper ground truth data
        [quality_metrics, aligned_results] = computeQualityMetrics(...
            Aout, Xout, bout, Y, ...
            A0_eval, X0_eval);  % Pass ground truth data directly

        % Store everything in extras
        extras.quality_metrics = quality_metrics;
        extras.aligned_results = aligned_results;
        extras.param_combinations = parameter_combinations;
        extras.param_idx = param_idx;
        extras.dataset_idx = dataset_idx;
        
    catch ME
        rethrow(ME);
    end
end