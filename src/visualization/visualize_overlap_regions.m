function visualize_overlap_regions(Xiter, k)
    % Get dimensions
    num_kernels = size(Xiter, 3);
    
    % Calculate overlap regions
    IDX = calculate_overlap_regions(Xiter, k);
    
    % Generate distinct colors for each kernel using lines or hsv
    colors = lines(num_kernels);
    
    % Visualize
    figure('Position', [100, 100, 300*num_kernels, 500]);
    
    % Plot activation regions (first row)
    for n = 1:num_kernels
        subplot(2, num_kernels, n);
        
        % Create kernel of ones
        A = ones(k(n, 1), k(n, 2));
        activation_region = conv2(Xiter(:,:,n), A, 'same');
        activation_region = activation_region / max(activation_region(:));
        
        imagesc(activation_region);
        title(sprintf('Kernel %d Activation Region', n));
        colorbar;
        axis equal tight;
    end
    
    % Plot overlap regions (second row)
    for n = 1:num_kernels
        subplot(2, num_kernels, n + num_kernels);
        
        % Create empty image for overlaps
        overlap_image = zeros(size(Xiter(:,:,1)));
        
        % Plot red dots where Xiter is non-zero
        [rows, cols] = find(Xiter(:,:,n));
        scatter(cols, rows, 'r.', 'SizeData', 50);
        axis ij;
        hold on;
        
        % Add overlaps with different colors for each kernel
        other_kernels = setdiff(1:num_kernels, n);  % Get indices of other kernels
        for idx = 1:length(other_kernels)
            m = other_kernels(idx);
            overlap_image(IDX(:,:,n,m)) = idx;  % Set discrete values
        end
        
        % Display overlap image
        h = imagesc(overlap_image);
        set(h, 'AlphaData', 0.3 * (overlap_image > 0));
        
        % Create custom colormap for this subplot
        other_colors = colors(other_kernels,:);  % Get colors for other kernels
        cmap = [1 1 1; other_colors];  % White background + colors for other kernels
        colormap(gca, cmap);
        
        % Create labels for colorbar including background
        tick_labels = cell(length(other_kernels) + 1, 1);
        tick_labels{1} = 'Background';
        for idx = 1:length(other_kernels)
            tick_labels{idx+1} = sprintf('K%d', other_kernels(idx));
        end
        
        title(sprintf('Kernel %d Overlaps', n));
        
        % Set up discrete colorbar
        c = colorbar;
        c.Ticks = 0:length(other_kernels);  % Include 0 for background
        c.TickLabels = tick_labels;
        c.Limits = [0, length(other_kernels)];
        
        axis equal tight;
        hold off;
    end
end