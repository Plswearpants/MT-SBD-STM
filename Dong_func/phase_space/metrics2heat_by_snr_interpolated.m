function metrics2heat_by_snr_interpolated(dataset_metrics, metric_type, mode, interp_factor, density_scale, density_label)
% Create multiple 2D heatmaps with fixed SNR values.
% Each subplot shows defect density vs axis-3 variable with linearly interpolated colors.
%
% Inputs:
%   dataset_metrics : struct from loadMetricDataset_new
%   metric_type     : 'kernel' | 'combined' | 'multiplied' (default 'kernel')
%   mode            : 1 -> axis-3 is Nobs, 2 -> axis-3 is side_length_ratio (default 1)
%   interp_factor   : grid upsampling factor for interpolation (default 4)
%   density_scale   : convert x-axis density by dividing with this value (default 1)
%   density_label   : x-axis label text (default 'defect density')

    if nargin < 2 || isempty(metric_type)
        metric_type = 'kernel';
    end
    if nargin < 3 || isempty(mode)
        mode = 1;
    end
    if nargin < 4 || isempty(interp_factor)
        interp_factor = 4;
    end
    if nargin < 5 || isempty(density_scale)
        density_scale = 1;
    end
    if nargin < 6 || isempty(density_label)
        density_label = 'defect density';
    end

    switch mode
        case 1
            axis3_values = dataset_metrics.Nobs_values;
            axis3_label = 'N_{obs}';
            axis3_tick_fmt = '%.0f';
            [metric_data, metric_name] = get_metric_data_mode1(dataset_metrics, metric_type);
        case 2
            if ~isfield(dataset_metrics, 'side_length_ratio_values') || isempty(dataset_metrics.side_length_ratio_values)
                error('mode=2 requires dataset_metrics.side_length_ratio_values.');
            end
            axis3_values = dataset_metrics.side_length_ratio_values;
            axis3_label = 'side length ratio';
            axis3_tick_fmt = '%.3f';
            [metric_data, metric_name] = get_metric_data_mode2(dataset_metrics, metric_type);
        otherwise
            error('mode must be 1 (Nobs) or 2 (side_length_ratio).');
    end

    snr_values = dataset_metrics.SNR_values(:)';
    theta_values = dataset_metrics.theta_cap_values(:)';
    theta_plot_values = theta_values ./ density_scale;

    num_snr = numel(snr_values);
    num_cols = ceil(sqrt(num_snr));
    num_rows = ceil(num_snr / num_cols);

    figure('Position', [100 100 420 * num_cols 420 * num_rows]);
    for s = 1:num_snr
        subplot(num_rows, num_cols, s);
        metric_slice = squeeze(metric_data(s, :, :)); % [theta x axis3]
        plot_interpolated_slice(theta_plot_values, axis3_values, metric_slice, ...
            sprintf('%s | SNR=%.2f', metric_name, snr_values(s)), axis3_label, axis3_tick_fmt, interp_factor, density_label);
    end

    sgtitle(sprintf('%s: interpolated 2D heatmaps by SNR', metric_name), ...
        'FontSize', 14, 'FontWeight', 'bold');
end

function plot_interpolated_slice(theta_values, axis3_values, metric_slice, title_str, axis3_label, axis3_tick_fmt, interp_factor, density_label)
    % Interpolate on a denser grid with linear interpolation.
    theta_log = log10(theta_values(:)');
    y_vals = axis3_values(:)';

    [X, Y] = meshgrid(theta_log, y_vals);     % base grid
    Z = metric_slice';                        % [axis3 x theta]

    xq = linspace(min(theta_log), max(theta_log), max(2, interp_factor * numel(theta_values)));
    yq = linspace(min(y_vals), max(y_vals), max(2, interp_factor * numel(y_vals)));
    [Xq, Yq] = meshgrid(xq, yq);
    Zq = interp2(X, Y, Z, Xq, Yq, 'linear');

    imagesc(10.^xq, yq, Zq);
    set(gca, 'YDir', 'normal');

    try
        cmap_data = load('metric_colormapv3.mat');
        colormap(cmap_data.CustomColormap);
    catch
        colormap('parula');
    end

    clim([0 1]);
    cb = colorbar;
    ylabel(cb, 'Performance Score');

    xlabel(density_label);
    ylabel(axis3_label);
    title(title_str, 'FontSize', 10);

    ax = gca;
    ax.XScale = 'log';
    ax.YScale = 'linear';
    ax.XTick = theta_values;
    ax.YTick = axis3_values;
    ax.XTickLabel = arrayfun(@(x) sprintf('%.0e', x), theta_values, 'UniformOutput', false);
    ax.YTickLabel = arrayfun(@(x) sprintf(axis3_tick_fmt, x), axis3_values, 'UniformOutput', false);
    ax.XTickLabelRotation = 45;
    axis square tight;
    grid on;
    ax.GridAlpha = 0.3;
end

function [metric_data, metric_name] = get_metric_data_mode1(dataset_metrics, metric_type)
    switch lower(metric_type)
        case 'kernel'
            metric_data = average_over_repetitions(dataset_metrics.kernel_quality_final);
            metric_name = 'Kernel Similarity';
        case 'combined'
            metric_data = average_over_repetitions(dataset_metrics.combined_activationScore);
            metric_name = 'Combined Activation Score';
        case 'multiplied'
            k = average_over_repetitions(dataset_metrics.kernel_quality_final);
            a = average_over_repetitions(dataset_metrics.activation_similarity_final);
            metric_data = k .* a;
            metric_name = 'Kernel x Activation';
        otherwise
            error('metric_type must be ''kernel'', ''combined'', or ''multiplied''.');
    end
end

function [metric_data, metric_name] = get_metric_data_mode2(dataset_metrics, metric_type)
    if ~isfield(dataset_metrics, 'kernel_quality_final_by_side_length_ratio')
        error('mode=2 requires *_by_side_length_ratio fields.');
    end

    switch lower(metric_type)
        case 'kernel'
            metric_data = average_over_repetitions(dataset_metrics.kernel_quality_final_by_side_length_ratio);
            metric_name = 'Kernel Similarity';
        case 'combined'
            metric_data = average_over_repetitions(dataset_metrics.combined_activationScore_by_side_length_ratio);
            metric_name = 'Combined Activation Score';
        case 'multiplied'
            k = average_over_repetitions(dataset_metrics.kernel_quality_final_by_side_length_ratio);
            a = average_over_repetitions(dataset_metrics.activation_similarity_final_by_side_length_ratio);
            metric_data = k .* a;
            metric_name = 'Kernel x Activation';
        otherwise
            error('metric_type must be ''kernel'', ''combined'', or ''multiplied''.');
    end
end

function avg_data = average_over_repetitions(data)
    if ndims(data) == 4
        avg_data = mean(data, 4, 'omitnan');
    else
        avg_data = data;
    end
end
