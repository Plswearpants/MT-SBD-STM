function metrics2heat_by_defect_density(dataset_metrics, metric_type, mode)
    % Create multiple 2D heatmap plots with fixed defect density values.
    % Each subplot shows SNR vs axis-3 variable with metric as color.
    %
    % Inputs:
    %   dataset_metrics : struct from loadMetricDataset_new
    %   metric_type     : 'kernel' | 'combined' | 'multiplied' (default 'kernel')
    %   mode            : 1 -> axis-3 is Nobs (default), 2 -> axis-3 is side_length_ratio
    
    if nargin < 2 || isempty(metric_type)
        metric_type = 'kernel';
    end
    if nargin < 3 || isempty(mode)
        mode = 1;
    end
    loader_mode = get_loader_axis_mode(dataset_metrics);
    
    switch mode
        case 1
            axis3_values = dataset_metrics.Nobs_values;
            axis3_label = 'N_{obs}';
            axis3_tick_fmt = '%.0f';
            [metric_data, metric_name] = get_metric_data(dataset_metrics, metric_type, false, loader_mode);
        case 2
            if ~isfield(dataset_metrics, 'side_length_ratio_values') || isempty(dataset_metrics.side_length_ratio_values)
                error('mode=2 requires dataset_metrics.side_length_ratio_values.');
            end
            axis3_values = dataset_metrics.side_length_ratio_values;
            axis3_label = 'side length ratio';
            axis3_tick_fmt = '%.3f';
            [metric_data, metric_name] = get_metric_data(dataset_metrics, metric_type, true, loader_mode);
        otherwise
            error('mode must be 1 (Nobs) or 2 (side_length_ratio).');
    end
    
    SNR_values = dataset_metrics.SNR_values;
    theta_cap_values = dataset_metrics.theta_cap_values;
    
    num_theta = length(theta_cap_values);
    num_cols = ceil(sqrt(num_theta));
    num_rows = ceil(num_theta / num_cols);
    
    figure('Position', [100 100 400*num_cols 400*num_rows]);
    
    for theta_idx = 1:num_theta
        theta_val = theta_cap_values(theta_idx);
        metric_slice = squeeze(metric_data(:, theta_idx, :));
        
        subplot(num_rows, num_cols, theta_idx);
        plot_2D_slice(SNR_values, axis3_values, metric_slice, ...
            sprintf('%s\nDefect Density: %.2e', metric_name, theta_val), ...
            axis3_label, axis3_tick_fmt);
    end
    
    sgtitle(sprintf('%s vs SNR and %s (Fixed Defect Density)', metric_name, axis3_label), ...
        'FontSize', 14, 'FontWeight', 'bold');
end

function [metric_data, metric_name] = get_metric_data(dataset_metrics, metric_type, use_side, loader_mode)
    suffix = '';
    if use_side && loader_mode ~= 2
        suffix = '_by_side_length_ratio';
    end
    
    switch lower(metric_type)
        case 'kernel'
            metric_data = average_over_repetitions(select_field(dataset_metrics, ...
                ['kernel_quality_final' suffix], 'kernel_quality_final'));
            metric_name = 'Kernel Similarity';
        case 'combined'
            metric_data = average_over_repetitions(select_field(dataset_metrics, ...
                ['combined_activationScore' suffix], 'combined_activationScore'));
            metric_name = 'Combined Activation Score';
        case 'multiplied'
            k = average_over_repetitions(select_field(dataset_metrics, ...
                ['kernel_quality_final' suffix], 'kernel_quality_final'));
            a = average_over_repetitions(select_field(dataset_metrics, ...
                ['activation_similarity_final' suffix], 'activation_similarity_final'));
            metric_data = k .* a;
            metric_name = 'Kernel \times Activation';
        otherwise
            error('metric_type must be ''kernel'', ''combined'', or ''multiplied''.');
    end
end

function data = select_field(metrics, preferred_name, fallback_name)
    if isfield(metrics, preferred_name)
        data = metrics.(preferred_name);
    elseif isfield(metrics, fallback_name)
        data = metrics.(fallback_name);
    else
        error('Missing metric field(s): %s and %s', preferred_name, fallback_name);
    end
end

function loader_mode = get_loader_axis_mode(dataset_metrics)
    if isfield(dataset_metrics, 'axis3_mode') && isscalar(dataset_metrics.axis3_mode)
        loader_mode = dataset_metrics.axis3_mode;
    else
        loader_mode = 1;
    end
end

function h = plot_2D_slice(SNR_values, axis3_values, metric_slice, title_str, axis3_label, axis3_tick_fmt)
    metric_display = metric_slice';
    
    h = imagesc(axis3_values, SNR_values, metric_display);
    set(gca, 'YDir', 'normal');
    
    try
        colormap(slanCM('viridis'));
    catch
        colormap('parula');
    end
    
    colorbar;
    clim([0 1]);
    
    xlabel(axis3_label);
    ylabel('SNR');
    title(title_str, 'FontSize', 10);
    
    cb = colorbar;
    ylabel(cb, 'Performance Score');
    
    ax = gca;
    ax.XScale = 'linear';
    ax.YScale = 'linear';
    axis square;
    ax.XTick = axis3_values;
    ax.YTick = SNR_values;
    
    ax.XTickLabel = arrayfun(@(x) sprintf(axis3_tick_fmt, x), axis3_values, 'UniformOutput', false);
    ax.YTickLabel = arrayfun(@(x) sprintf('%.1f', x), SNR_values, 'UniformOutput', false);
    ax.XTickLabelRotation = 45;
    
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

