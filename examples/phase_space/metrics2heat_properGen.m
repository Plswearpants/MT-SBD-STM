function metrics2heat_properGen(dataset_metrics)
    % Create 3D heatspace visualization for both kernel similarity and combined activation metrics
    % Input:
    %   dataset_metrics: structure containing:
    %       - SNR_values
    %       - theta_cap_values
    %       - Nobs_values
    %       - kernel_similarity [SNR × theta_cap × Nobs]
    %       - activation_metrics [SNR × theta_cap × Nobs]
    
    % Create figure with two subplots
    figure('Position', [100 100 1600 800]);
    
    % 1. Kernel Similarity Heatspace
    subplot(1,2,1);
    h1 = plot_3D_heatspace(dataset_metrics.SNR_values, ...
                          dataset_metrics.theta_cap_values, ...
                          dataset_metrics.Nobs_values, ...
                          dataset_metrics.kernel_quality_final, ...
                          'Kernel Similarity');
    
    % 2. Combined Activation Score Heatspace
    subplot(1,2,2);
    h2 = plot_3D_heatspace(dataset_metrics.SNR_values, ...
                          dataset_metrics.theta_cap_values, ...
                          dataset_metrics.Nobs_values, ...
                          dataset_metrics.combined_activationScore, ...
                          'Combined Activation Score');
end

function h = plot_3D_heatspace(SNR_values, theta_cap_values, Nobs_values, metric_values, title_str)
    % Create meshgrid matching the data structure: Metric(SNR, theta_cap, Nobs)
    [SNR, Theta, Nobs] = meshgrid(SNR_values, theta_cap_values, Nobs_values);
    
    % Ensure all inputs are column vectors
    x = Theta(:);   % Column vector for x coordinates (theta_cap)
    y = Nobs(:);    % Column vector for y coordinates (Nobs)
    z = SNR(:);     % Column vector for z coordinates (SNR)
    c = metric_values(:);  % Column vector for colors

    % Create scatter plot with color based on metric values
    h = scatter3(x, y, z, 100, c, 'filled');
    
    % Load custom colormap
    cmap_data = load('custom_PiYG_colormap.mat');
    colormap(cmap_data.custom_map);
    
    % Customize appearance
    colorbar;
    clim([0 1]);
    
    % Add labels and title
    xlabel('defect density');
    ylabel('N_{obs}');
    zlabel('SNR');
    title(title_str);
    
    % Add grid and adjust view
    grid on;
    view(-35, 78);
    
    % Add colorbar label
    c = colorbar;
    ylabel(c, 'Performance Score');
    
    % Make axis labels more readable and set directions
    ax = gca;
    ax.XScale = 'log';  % Theta cap in log scale
    ax.YScale = 'linear';  % Nobs in linear scale
    ax.ZScale = 'linear';  % SNR in linear scale
    
    % Set axis directions
    set(ax, 'XDir', 'normal');   % theta_cap: low -> high
    set(ax, 'YDir', 'reverse');  % Nobs: high -> low
    set(ax, 'ZDir', 'reverse');  % SNR: high -> low
    
    % Set tick values to match the actual data points
    ax.XTick = theta_cap_values;
    ax.YTick = Nobs_values;
    ax.ZTick = SNR_values;
    
    % Format tick labels
    ax.XTickLabel = arrayfun(@(x) sprintf('%.0e', x), theta_cap_values, 'UniformOutput', false);
    ax.YTickLabel = arrayfun(@(x) sprintf('%.0f', x), Nobs_values, 'UniformOutput', false);
    ax.ZTickLabel = arrayfun(@(x) sprintf('%.1f', x), SNR_values, 'UniformOutput', false);
    
    % Rotate x-axis labels for better readability
    ax.XTickLabelRotation = 45;
end
