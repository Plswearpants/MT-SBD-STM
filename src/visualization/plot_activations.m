function plot_activations(activation_array)
    % PLOT_ACTIVATIONS Plot activations with different colors for each n value
    %   activation_array: 3D array of shape (x,y,n) where n represents different activation types
    
    % Get dimensions
    [height, width, num_types] = size(activation_array);
    
    % Create figure
    figure('Name', 'Activation Plot', 'Position', [100 100 800 600]);
    
    % Generate distinct colors for each n value
    colors = lines(num_types);
    
    % Create legend entries
    legend_entries = cell(num_types, 1);
    
    % Plot each activation type with a different color
    hold on;
    for n = 1:num_types
        % Get current activation map
        current_activation = activation_array(:,:,n);
        
        % Calculate 80th percentile threshold
        threshold = prctile(current_activation(current_activation > 0.01), 80);

        % Find activations above threshold
        [y, x] = find(current_activation > threshold);
        
        if ~isempty(x)
            % Plot points for this n value
            scatter(x, y, 50, colors(n,:), 'filled', 'o');
            
            % Add to legend with threshold value
            legend_entries{n} = sprintf('n = %d (threshold = %.2f)', n, threshold);
        end
    end
    hold off;
    
    % Add labels and title
    xlabel('X');
    ylabel('Y');
    title('Activation Plot (Points above 80th percentile)');
    
    % Add legend (only for non-empty activation types)
    legend(legend_entries(~cellfun(@isempty, legend_entries)), 'Location', 'best');
    
    % Set equal aspect ratio and adjust axis limits
    axis equal;
    xlim([1 width]);
    ylim([1 height]);
    
    % Add grid
    grid on;
end 