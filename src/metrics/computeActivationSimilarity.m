function [similarity, filtered_maps] = computeActivationSimilarity(X0, Xout, kernel_size, visualize, indices)
    % Compute similarity between original and reconstructed activations using 
    % density-adaptive Gaussian filtering
    %
    % Inputs:
    %   X0: original activation maps [height x width x num_kernels]
    %   Xout: reconstructed activation maps [height x width x num_kernels]
    %   kernel_size: size of each kernel [num_kernels x 2]
    %   visualize: (optional) boolean for visualization
    %   indices: Optional [dataset_num, param_num] for labeling plots
    %
    % Outputs:
    %   similarity: array of similarity scores for each kernel
    %   filtered_maps: struct array containing filtered maps and parameters
    
    if nargin < 4
        visualize = false;
    end
    
    % Prepare index string if indices provided
    idx_str = '';
    if nargin > 4 && ~isempty(indices)
        idx_str = sprintf('\nDataset # %d, Parameter Set # %d', indices(1), indices(2));
    end
    
    num_kernels = size(X0, 3);
    similarity = zeros(1, num_kernels);
    
    % Process each kernel
    for k = 1:num_kernels
        % Compute density-adaptive sigma
        [sigma, L, density] = computeAdaptiveSigma(X0(:,:,k), kernel_size(k,:));
        
        % Create Gaussian filter
        window_size = ceil(3 * sigma);
        [x, y] = meshgrid(-window_size:window_size);
        gaussian_kernel = exp(-(x.^2 + y.^2)/(2*sigma^2));
        gaussian_kernel = gaussian_kernel / sum(gaussian_kernel(:));
        
        % Apply Gaussian filter
        X0_filtered = conv2(X0(:,:,k), gaussian_kernel, 'same');
        Xout_filtered = conv2(Xout(:,:,k), gaussian_kernel, 'same');
        
        % Compute normalized dot product
        similarity(k) = dot(X0_filtered(:), Xout_filtered(:)) / ...
                       (norm(X0_filtered(:)) * norm(Xout_filtered(:)));
        
        % Store results for this kernel
        filtered_maps(k).X0 = X0_filtered;
        filtered_maps(k).Xout = Xout_filtered;
        filtered_maps(k).gaussian_kernel = gaussian_kernel;
        filtered_maps(k).sigma = sigma;
        filtered_maps(k).L = L;
        filtered_maps(k).density = density;
        
        % Visualize if requested
        if visualize       
            visualizeSimilarityAnalysis(X0(:,:,k), Xout(:,:,k), ...
                X0_filtered, Xout_filtered, gaussian_kernel, ...
                similarity(k), sigma, L, density, k, num_kernels, visualize, idx_str);
        end
    end
    
    % Plot summary if visualizing
    %   if visualize && num_kernels > 1
    %    plotSimilaritySummary(similarity, filtered_maps, idx_str);
    %   end
end

function [sigma, L, density] = computeAdaptiveSigma(X0, kernel_size)
    % Compute adaptive sigma based on both kernel size and activation density
    %
    % For a 2D Poisson process with density ρ:
    % - L is chosen such that P(no other activations within radius L) = 0.95
    % - This gives: L = sqrt(-ln(0.95)/(ρπ))
    % - σ is then chosen as min(L/3, kernel_size/10)
    
    density = sum(X0(:) > 0) / numel(X0);
    L = sqrt(-log(0.95)/(density * pi));
    
    sigma_kernel = min(kernel_size)/10;
    sigma_density = L/3;
    sigma = min(sigma_kernel, sigma_density);
end

function visualizeSimilarityAnalysis(X0, Xout, X0_filtered, Xout_filtered, ...
                                   gaussian_kernel, similarity, sigma, L, density, ...
                                   k, num_kernels, visualize, idx_str)
    % Create figure for first kernel or if only one kernel if visualize is true
    if visualize
        if k == 1 || num_kernels == 1
            figure('Name', ['Activation Similarity Analysis' idx_str], ...
               'Position', [100, 100, 1200, 300*num_kernels]);
        end
        
        % Calculate subplot positions
        base_idx = (k-1)*4;
        % Gaussian kernel
        subplot(num_kernels,4,base_idx + 1);
        surf(gaussian_kernel);
        title(sprintf('K%d: Filter\nσ=%.1f, L=%.1f', k, sigma, L));
        axis square;
        
        % Filtered maps
        subplot(num_kernels,4,base_idx + 2);
        imagesc(X0_filtered);
        title(sprintf('K%d: Filtered X0', k));
        colorbar; axis square;
        
        subplot(num_kernels,4,base_idx + 3);
        imagesc(Xout_filtered);
        title(sprintf('K%d: Filtered Xout', k));
        colorbar; axis square;
        
        % Difference
        subplot(num_kernels,4,base_idx + 4);
        imagesc(X0_filtered - Xout_filtered);
        title(sprintf('K%d: Diff\nSim=%.3f', k, similarity));
        colorbar; axis square;

        sgtitle(['Kernel Similarity Analysis' idx_str]);
    end
end

function plotSimilaritySummary(similarity, filtered_maps, idx_str)
    figure('Name', ['Similarity Summary' idx_str]);
    
    % Plot similarity scores
    subplot(2,2,1);
    bar(similarity);
    title('Similarity Scores');
    xlabel('Kernel'); ylabel('Score');
    ylim([0 1]);
    
    % Plot densities
    subplot(2,2,2);
    densities = [filtered_maps.density];
    bar(densities);
    title('Activation Densities');
    xlabel('Kernel'); ylabel('Density');
    set(gca, 'YScale', 'log');
    
    % Plot L values
    subplot(2,2,3);
    L_values = [filtered_maps.L];
    bar(L_values);
    title('Isolation Lengths (L)');
    xlabel('Kernel'); ylabel('Pixels');
    
    % Plot sigma values
    subplot(2,2,4);
    sigma_values = [filtered_maps.sigma];
    bar(sigma_values);
    title('Filter Widths (σ)');
    xlabel('Kernel'); ylabel('Pixels');
    
    sgtitle(['Summary Statistics Across Kernels' idx_str]);
end
