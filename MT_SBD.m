function [ Aout, Xout, bout, extras ] = MT_SBD( Y, k, params, dispfun, kernel_initialguess, Max_iteration, maxIT)
    %Multi Type Sparse Blind Deconvolution for real data 
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
    
    %% Process input arguments
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
    
    if ~isfield(params, 'xinit') || isempty(params.xinit)
        xinit = [];
    else
        xinit = params.xinit;
    end
    
    if ~isfield(params, 'noise_var')
        error('params must contain noise_var for quality metrics');
    end
    noise_var = params.noise_var;
    
    kernel_num = size(k,1);
    mu = 10^-6;
    
    % Update the configuration file with the new max_iteration
    update_config('Xsolve_config.mat', 'MAXIT', Max_iteration, 'Xsolve_config_tunable.mat');
    update_config('Asolve_config.mat','options.maxiter', Max_iteration, 'Asolve_config_tunable.mat');
    
    %% Phase I: Initialization and First Iteration with demixing applied
    fprintf('PHASE I: Initialization and First Iteration\n');
    A = kernel_initialguess;
    X_struct = struct();
    Xiter = zeros([size(Y),kernel_num]);
    biter = zeros(kernel_num,1);
    
    % Initialize new metrics arrays instead of the old ones
    extras.phase1.quality_metrics = zeros(maxIT, 1);
    extras.phase1.residuals = zeros([size(Y), maxIT]);
    
    % Add demixing factor
    faint_factor = 1;
    
    % Compute initial kernel order based on variance of initial guess
    kernel_variances = zeros(1, kernel_num);
    for n = 1:kernel_num
        kernel_variances(n) = var(kernel_initialguess{n}(:));
    end
    [~, kernel_order] = sort(kernel_variances, 'descend');
    
    % Print the initial kernel processing order
    fprintf('Initial kernel processing order: ');
    fprintf('%d ', kernel_order);
    fprintf('\n');
    
    % Print the initial variance of each kernel
    fprintf('Initial kernel variances: ');
    for n = 1:kernel_num
        fprintf('%.6f ', kernel_variances(n));
    end
    fprintf('\n');
    
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
        
        % Update each kernel in the current order
        for idx = 1:kernel_num
            n = kernel_order(idx);
            % Calculate Yiter for this kernel (changed to demixing approach)
            Yiter = Y_residual + (iter > 1) * (1-1/(faint_factor*iter+1))*convfft2(A{n}, Xiter(:,:,n)) + (1/(faint_factor*iter+1))*Y_sum;
            Y_residual_pre = Y_residual + convfft2(A{n}, Xiter(:,:,n));
            dispfun1 = @(A, X) dispfun{n}(Y, A, X, k(n,:), []);
            
            if iter == 1
                if isempty(xinit)
                    % Initial X computation without xinit
                    X_struct.(['x',num2str(n)]) = Xsolve_FISTA_tunable(Y, A{n}, lambda1(n), mu, [], xpos);
                else
                    % Initial X computation with xinit
                    X_struct.(['x',num2str(n)]) = Xsolve_FISTA_tunable(Y, A{n}, lambda1(n), mu, xinit{n}, xpos);
                end
            end
            
            [A{n}, X_struct.(['x',num2str(n)]), info] = Asolve_Manopt_tunable(Yiter, A{n}, lambda1(n), Xsolve, X_struct.(['x',num2str(n)]), xpos, getbias, dispfun1);
            
            Xiter(:,:,n) = X_struct.(['x',num2str(n)]).X;
            biter(n) = X_struct.(['x',num2str(n)]).b;
            Y_residual = Y_residual_pre - convfft2(A{n},Xiter(:,:,n));
        end
    
        % Compute quality metrics using residual
        [quality_metric, residual] = computeResidualQuality(Y, A, Xiter, noise_var);
        extras.phase1.quality_metrics(iter) = quality_metric;
        extras.phase1.residuals(:,:,iter) = residual;
        
        iter_runtime = toc(iter_starttime);
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
    
        lambda = lambda1;
        lam2fac = (lambda2./lambda1).^(1/nrefine);
        
        % Initialize Phase II metrics
        extras.phase2.quality_metrics = zeros(nrefine + 1, 1);
        extras.phase2.residuals = zeros([size(Y), nrefine + 1]);
        
        for i = 1:nrefine + 1
            fprintf('lambda iteration %d/%d: \n', i, nrefine + 1);
            
            % Initialize Y_sum to zero before accumulation
            Y_sum = zeros(size(Y));
            % Standard Y_residual calculation (no demixing)
            for m = 1:kernel_num
                Y_sum = Y_sum + convfft2(A2{m}, X2_struct.(['x',num2str(m)]).X);
            end
            Y_residual = Y - Y_sum;
    
            for idx = 1:kernel_num
                n = kernel_order(idx);
                fprintf('Processing kernel %d, lambda = %.1e: \n', n, lambda(n));
                % Calculate Yiter without demixing
                Yiter = Y_residual + convfft2(A2{n}, X2_struct.(['x',num2str(n)]).X);
                Y_residual_pre = Y_residual + convfft2(A2{n}, X2_struct.(['x',num2str(n)]).X);
                
                dispfun2 = @(A, X) dispfun{n}(Y, A, X, k3(n,:), kplus(n,:));
                [A2{n}, X2_struct.(['x',num2str(n)]), info] = Asolve_Manopt_tunable(Yiter, A2{n}, lambda(n), Xsolve, X2_struct.(['x',num2str(n)]), xpos, getbias, dispfun2);
                
                % Attempt to 'unshift' the a and x
                score = zeros(2*kplus(n,1)+1, 2*kplus(n,2)+1);
                % Calculate center position
                center_y = kplus(n,1) + 1;
                center_x = kplus(n,2) + 1;
                
                for tau1 = -kplus(n,1):kplus(n,1)
                    ind1 = tau1+kplus(n,1)+1;
                    for tau2 = -kplus(n,2):kplus(n,2)
                        ind2 = tau2+kplus(n,2)+1;
                        
                        % Get the selected region
                        temp = A2{n}(ind1:(ind1+k(n,1)-1), ind2:(ind2+k(n,2)-1));
                        
                        % Calculate decay for each point in temp relative to its center
                        [y_coords, x_coords] = ndgrid(1:size(temp,1), 1:size(temp,2));
                        center_y_temp = (size(temp,1)+1)/2;
                        center_x_temp = (size(temp,2)+1)/2;
                        
                        % Calculate normalized distances from center (0 to 1)
                        dist_y = abs(y_coords - center_y_temp) / center_y_temp;
                        dist_x = abs(x_coords - center_x_temp) / center_x_temp;
                        dist = sqrt(dist_y.^2 + dist_x.^2);
                        
                        % Calculate decay factor (1 at center, 0.5 at edges)
                        decay = exp(-log(2) * dist);  % exp(-log(2)) = 0.5 at dist=1
                        
                        % Apply decay to temp
                        temp_decayed = temp .* decay;
                        score(ind1,ind2) = norm(temp_decayed(:), 1);
                    end
                end
                [temp,ind1] = max(score); [~,ind2] = max(temp);
                tau = [ind1(ind2) ind2]-kplus(n,:)-1;
                A2{n} = circshift(A2{n},-tau);
                X2_struct.(['x',num2str(n)]).X = circshift(X2_struct.(['x',num2str(n)]).X,tau);
                X2_struct.(['x',num2str(n)]).W = circshift(X2_struct.(['x',num2str(n)]).W,tau);
    
                % Update Y_residual after processing this kernel
                Y_residual = Y_residual_pre - convfft2(A2{n}, X2_struct.(['x',num2str(n)]).X);
    
                % Save phase 2 extras
                extras.phase2.A{n} = A2{n};
                extras.phase2.X{n} = X2_struct.(['x',num2str(n)]).X;
                extras.phase2.b(n) = X2_struct.(['x',num2str(n)]).b;
                extras.phase2.info{n} = info;
            end
            
            % Compute quality metrics using residual
            X2_combined = zeros(size(Y,1), size(Y,2), kernel_num);
            A2_central = cell(1, kernel_num);
            for n = 1:kernel_num
                X2_combined(:,:,n) = X2_struct.(['x',num2str(n)]).X;
                A2_central{n} = A2{n}(kplus(n,1)+(1:k(n,1)), kplus(n,2)+(1:k(n,2)));
            end
    
            [quality_metric, residual] = computeResidualQuality(Y, A2_central, X2_combined, noise_var);
            extras.phase2.quality_metrics(i) = quality_metric;
            extras.phase2.residuals(:,:,i) = residual;
            
            % Update the results printing section
            fprintf('Refinement %d Metrics:\n', i);
            fprintf('Quality Metric: %.3f\n', quality_metric);
            
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
    
    % Store runtime in extras
    extras.runtime = total_runtime;
end