function metrics2heat_by_defect_density(dataset_metrics, metric_type)
    % Create multiple 3D plots with fixed defect density values
    % Each plot shows SNR vs Nobs with metric as color
    %
    % Input:
    %   dataset_metrics: structure containing:
    %       - SNR_values
    %       - theta_cap_values (defect density values)
    %       - Nobs_values
    %       - kernel_quality_final [SNR × theta_cap × Nobs] or [SNR × theta_cap × Nobs × rep]
    %       - combined_activationScore [SNR × theta_cap × Nobs] or [SNR × theta_cap × Nobs × rep]
    %   metric_type: 'kernel' or 'combined' (default: 'kernel')
    %                'kernel' -> kernel_quality_final
    %                'combined' -> combined_activationScore
    %
    % Output: Creates figure(s) with subplots, one for each defect density value
    
    if nargin < 2
        metric_type = 'kernel';
    end
    
    % Average over repetitions if 4D arrays are present
    if strcmpi(metric_type, 'kernel')
        metric_data = average_over_repetitions(dataset_metrics.kernel_quality_final);
        metric_name = 'Kernel Similarity';
    elseif strcmpi(metric_type, 'combined')
        metric_data = average_over_repetitions(dataset_metrics.combined_activationScore);
        metric_name = 'Combined Activation Score';
    else
        error('metric_type must be ''kernel'' or ''combined''');
    end
    
    % Get parameter values
    SNR_values = dataset_metrics.SNR_values;
    theta_cap_values = dataset_metrics.theta_cap_values;
    Nobs_values = dataset_metrics.Nobs_values;
    
    num_theta = length(theta_cap_values);
    
    % Determine subplot layout (try to make it roughly square)
    num_cols = ceil(sqrt(num_theta));
    num_rows = ceil(num_theta / num_cols);
    
    % Create figure
    fig = figure('Position', [100 100 400*num_cols 400*num_rows]);
    
    % Loop over each defect density value
    for theta_idx = 1:num_theta
        theta_val = theta_cap_values(theta_idx);
        
        % Extract 2D slice for this defect density: [SNR × Nobs]
        metric_slice = squeeze(metric_data(:, theta_idx, :));
        
        % Create subplot
        subplot(num_rows, num_cols, theta_idx);
        
        % Create 3D plot for this defect density
        plot_3D_slice(SNR_values, Nobs_values, metric_slice, ...
            sprintf('%s\nDefect Density: %.2e', metric_name, theta_val));
    end
    
    % Add overall title
    sgtitle(sprintf('%s vs SNR and N_{obs} (Fixed Defect Density)', metric_name), ...
        'FontSize', 14, 'FontWeight', 'bold');
end

function h = plot_3D_slice(SNR_values, Nobs_values, metric_slice, title_str)
    % Create 2D heatmap plot for a 2D metric slice
    % Input:
    %   SNR_values: vector of SNR values
    %   Nobs_values: vector of Nobs values
    %   metric_slice: 2D array [SNR × Nobs] with metric values
    %   title_str: title for the plot
    
    % Use imagesc or pcolor for 2D heatmap
    % Note: imagesc expects [rows × cols] where rows correspond to Y-axis
    % Our data is [SNR × Nobs], so we need to transpose for display
    % imagesc displays with first dimension as Y (rows) and second as X (cols)
    
    % Transpose so that SNR is on Y-axis and Nobs is on X-axis
    metric_display = metric_slice';
    
    % Create the heatmap
    h = imagesc(Nobs_values, SNR_values, metric_display);
    set(gca, 'YDir', 'normal');  % Normal Y direction (low to high)
    
    % Load custom colormap
    try
        cmap_data = load('metric_colormapv3.mat');
        colormap(cmap_data.CustomColormap);
    catch
        % Fallback to default colormap if file not found
        colormap('parula');
    end
    
    % Customize appearance
    colorbar;
    clim([0 1]);
    
    % Add labels and title
    xlabel('N_{obs}');
    ylabel('SNR');
    title(title_str, 'FontSize', 10);
    
    % Add colorbar label
    cb = colorbar;
    ylabel(cb, 'Performance Score');
    
    % Make axis labels more readable
    ax = gca;
    ax.XScale = 'linear';  % Nobs in linear scale
    ax.YScale = 'linear';  % SNR in linear scale
    
    % Set tick values
    ax.XTick = Nobs_values;
    ax.YTick = SNR_values;
    
    % Format tick labels
    ax.XTickLabel = arrayfun(@(x) sprintf('%.0f', x), Nobs_values, 'UniformOutput', false);
    ax.YTickLabel = arrayfun(@(x) sprintf('%.1f', x), SNR_values, 'UniformOutput', false);
    
    % Rotate x-axis labels for better readability
    ax.XTickLabelRotation = 45;
    
    % Add grid (optional, can be removed if too cluttered)
    grid on;
    ax.GridAlpha = 0.3;
end

function avg_data = average_over_repetitions(data)
    % Average over the 4th dimension (repetitions) if present
    % Input: data can be 3D [SNR × theta × N_obs] or 4D [SNR × theta × N_obs × rep]
    % Output: 3D array [SNR × theta × N_obs] averaged over repetitions
    
    if ndims(data) == 4
        % Average over the 4th dimension, ignoring NaN values
        avg_data = mean(data, 4, 'omitnan');
    elseif ndims(data) == 3
        % Already 3D, return as is
        avg_data = data;
    else
        % Unexpected dimensions, return as is
        avg_data = data;
    end
end

