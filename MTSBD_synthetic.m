function [ Aout, Xout, bout, extras ] = MTSBD_synthetic( Y, k, params, dispfun, kernel_initialguess, Max_iteration, maxIT)
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
    
    kernel_num = size(k,1);
    mu = 10^-6;
    
    % Update the configuration file with the new max_iteration
    update_config('Xsolve_config.mat', 'MAXIT', Max_iteration, 'Xsolve_config_tunable.mat');
    update_config('Asolve_config.mat','options.maxiter', Max_iteration, 'Asolve_config_tunable.mat');
    
    % Extract X0 and A0 from params
    if ~isfield(params, 'X0') || ~isfield(params, 'A0')
        error('params must contain X0 and A0 for quality metrics');
    end
    
    X0 = params.X0;
    A0 = params.A0;
    
    %% Phase I: Initialization and First Iteration with demixing applied
    fprintf('PHASE I: Initialization and First Iteration\n');
    A = kernel_initialguess;
    X_struct = struct();
    Xiter = zeros([size(Y),kernel_num]);
    biter = zeros(kernel_num,1);
    
    % Initialize quality metrics arrays in extras
    extras.phase1.activation_metrics = zeros(maxIT, kernel_num);
    extras.phase1.kernel_quality_factors = zeros(maxIT, kernel_num);
    
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
                % Initial X computation
                X_struct.(['x',num2str(n)]) = Xsolve_FISTA_tunable(Y, A{n}, lambda1(n), mu, xinit, xpos);
                else
                X_struct.(['x',num2str(n)]) = Xsolve_FISTA_tunable(Y, A{n}, lambda1(n), mu, xinit{n}, xpos);
                % Or initial A computation 
                end
           end
            
            [A{n}, X_struct.(['x',num2str(n)]), info] = Asolve_Manopt_tunable(Yiter, A{n}, lambda1(n), Xsolve, X_struct.(['x',num2str(n)]), xpos, getbias, dispfun1);
            
            Xiter(:,:,n) = X_struct.(['x',num2str(n)]).X;
            biter(n) = X_struct.(['x',num2str(n)]).b;
            Y_residual = Y_residual_pre - convfft2(A{n},Xiter(:,:,n));
        end
    
        % Keep the same quality metrics computation
        [activation_similarity, kernel_similarity] = computeQualityMetrics(X0, Xiter, A0, A, k);
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
        extras.phase2.activation_metrics = zeros(nrefine + 1, kernel_num);
        extras.phase2.kernel_quality_factors = zeros(nrefine + 1, kernel_num);
        extras.phase2.A = cell(1, kernel_num);
        extras.phase2.X = cell(1, kernel_num);
        extras.phase2.b = zeros(kernel_num, 1);
        extras.phase2.info = cell(1, kernel_num);
        
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
                
                % Attempt to 'unshift" the a and x
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
                        score(ind1,ind2) = norm(temp(:), 1);
                    end
                end
                [temp,ind1] = max(score); [~,ind2] = max(temp);
                tau = [ind1(ind2) ind2]-kplus(n,:)-1;
                
                % Visualize the process for the first iteration of Phase II
                if i == 1
                    % First figure for kernel visualization
                    figure('Position', [100, 100, 1500, 400]);
                    
                    % Plot original kernel
                    subplot(1,5,1);
                    imagesc(A{n});
                    %title(sprintf('Original Kernel %d', n));
                    colorbar;
                    axis equal tight;
                    
                    % Initialize padded kernel centered
                    k3 = k + 2*kplus;
                    A2_padded = zeros(k3);
                    % Place kernel at center
                    center_y = kplus(n,1) + 1;
                    center_x = kplus(n,2) + 1;
                    A2_padded(center_y+(1:k(n,1)), center_x+(1:k(n,2))) = A{n};
                    
                    % Plot padded kernel
                    subplot(1,5,2);
                    imagesc(A2_padded);
                    %title(sprintf('Padded Kernel %d', n));
                    colorbar;
                    axis equal tight;
                    
                    % Plot updated kernel before shifting
                    subplot(1,5,3);
                    imagesc(A2{n});
                    %title(sprintf('Updated Kernel %d (Asolve)', n));
                    colorbar;
                    axis equal tight;
                    
                    % Plot score heatmap
                    subplot(1,5,4);
                    imagesc(score);
                    %title(sprintf('Score Heatmap Kernel %d', n));
                    colorbar;
                    axis equal tight;
                    hold on;
                    plot(ind2, ind1(ind2), 'rx', 'MarkerSize', 10, 'LineWidth', 2);
                    hold off;
                    
                    % Plot final truncated kernel
                    A2_shifted = circshift(A2{n},-tau);
                    A2_truncated = A2_shifted(kplus(n,1)+(1:k(n,1)), kplus(n,2)+(1:k(n,2)));
                    subplot(1,5,5);
                    imagesc(A2_truncated);
                    %title(sprintf('Final Truncated Kernel %d', n));
                    colorbar;
                    axis equal tight;
                    
                    % Add overall title
                    %sgtitle(sprintf('Kernel %d Shifting Process (Phase II, Iteration 1)', n));
                    
                    % Second figure for activation map visualization
                    figure('Position', [100, 550, 1500, 400]);
                    
                    % Plot original activation map
                    subplot(1,4,1);
                    imagesc(X2_struct.(['x',num2str(n)]).X);
                    title(sprintf('Original Activation Map %d', n));
                    colorbar;
                    axis equal tight;
                    
                    % Plot activation map after shifting
                    X_shifted = circshift(X2_struct.(['x',num2str(n)]).X, tau);
                    subplot(1,4,2);
                    imagesc(X_shifted);
                    title(sprintf('Shifted Activation Map %d', n));
                    colorbar;
                    axis equal tight;
                    
                    % Plot Y_residual before and after
                    subplot(1,4,3);
                    imagesc(Y_residual_pre);
                    title(sprintf('Y Residual Before Shift %d', n));
                    colorbar;
                    axis equal tight;
                    
                    subplot(1,4,4);
                    imagesc(Y_residual);
                    title(sprintf('Y Residual After Shift %d', n));
                    colorbar;
                    axis equal tight;
                    
                    % Add overall title
                    sgtitle(sprintf('Activation Map and Residual Shifting Process (Phase II, Iteration 1)'));
                end
                
                % Apply shift
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
    
    % Store runtime in extras
    extras.runtime = total_runtime;
end