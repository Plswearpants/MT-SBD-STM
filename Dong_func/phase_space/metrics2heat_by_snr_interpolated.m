function metrics2heat_by_snr_interpolated(dataset_metrics, metric_type, mode, interp_factor, manual_colormap, contour_levels)
% Create multiple 2D heatmaps with fixed SNR values.
% Each subplot shows defect density vs axis-3 variable with linearly interpolated colors.
%
% Inputs:
%   dataset_metrics : struct from loadMetricDataset_new
%   metric_type     : 'kernel' | 'combined' | 'fidelity' | 'multiplied' | 'kernel_baseline' (default 'kernel')
%   mode            : 1 -> axis-3 is Nobs, 2 -> axis-3 is side_length_ratio (default 1)
%   interp_factor   : grid upsampling factor for interpolation (default 4)
%   manual_colormap : optional colormap override (e.g., slanCM('viridis'))
%   contour_levels  : optional vector of iso-metric levels (e.g., [0.85, 0.95])
%
    if nargin < 2 || isempty(metric_type)
        metric_type = 'kernel';
    end
    if nargin < 3 || isempty(mode)
        mode = 1;
    end
    if nargin < 4 || isempty(interp_factor)
        interp_factor = 4;
    end
    if nargin < 5
        manual_colormap = [];
    end
    if nargin < 6
        contour_levels = [];
    end
    use_unit_clim = ~strcmpi(metric_type, 'fidelity');
    loader_mode = get_loader_axis_mode(dataset_metrics);

    switch mode
        case 1
            axis3_values = dataset_metrics.Nobs_values;
            axis3_label = 'N_{obs}';
            axis3_tick_fmt = '%.0f';
            [metric_data, metric_name] = get_metric_data_mode1(dataset_metrics, metric_type, loader_mode);
        case 2
            if ~isfield(dataset_metrics, 'side_length_ratio_values') || isempty(dataset_metrics.side_length_ratio_values)
                error('mode=2 requires dataset_metrics.side_length_ratio_values.');
            end
            axis3_values = dataset_metrics.side_length_ratio_values;
            axis3_label = 'side length ratio';
            axis3_tick_fmt = '%.3f';
            [metric_data, metric_name] = get_metric_data_mode2(dataset_metrics, metric_type, loader_mode);
        otherwise
            error('mode must be 1 (Nobs) or 2 (side_length_ratio).');
    end

    snr_values = dataset_metrics.SNR_values(:)';
    theta_values = dataset_metrics.theta_cap_values(:)';
    metric_data = align_metric_tensor_layout(metric_data, numel(snr_values), numel(theta_values), numel(axis3_values), metric_type);

    num_snr = numel(snr_values);
    num_cols = ceil(sqrt(num_snr));
    num_rows = ceil(num_snr / num_cols);

    figure('Position', [100 100 420 * num_cols 420 * num_rows]);
    for s = 1:num_snr
        subplot(num_rows, num_cols, s);
        metric_slice = squeeze(metric_data(s, :, :)); % [theta x axis3]
        plot_interpolated_slice(theta_values, axis3_values, metric_slice, ...
            sprintf('%s | SNR=%.2f', metric_name, snr_values(s)), axis3_label, axis3_tick_fmt, ...
            interp_factor, manual_colormap, contour_levels, use_unit_clim);
    end

    sgtitle(sprintf('%s: interpolated 2D heatmaps by SNR', metric_name), ...
        'FontSize', 14, 'FontWeight', 'bold');
end

function plot_interpolated_slice(theta_values, axis3_values, metric_slice, title_str, axis3_label, axis3_tick_fmt, interp_factor, manual_colormap, contour_levels, use_unit_clim)
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

    apply_colormap(manual_colormap);

    if use_unit_clim
        clim([0 1]);
    else
        z_valid = Zq(~isnan(Zq));
        if ~isempty(z_valid)
            zmin = min(z_valid);
            zmax = max(z_valid);
            if zmin < zmax
                clim([zmin zmax]);
            end
        end
    end
    cb = colorbar;
    ylabel(cb, 'Performance Score');

    if ~isempty(contour_levels)
        contour_levels_plot = normalize_contour_levels(contour_levels);
        hold on;
        contour(10.^xq, yq, Zq, contour_levels_plot, ...
            'LineColor', 'k', 'LineStyle', '--', 'LineWidth', 1.2);
        hold off;
    end

    xlabel('defect density (per site)');
    ylabel(axis3_label);
    title(title_str, 'FontSize', 10);

    ax = gca;
    ax.YScale = 'linear';
    ax.YTick = axis3_values;
    apply_defect_density_tick_style(ax, theta_values);
    ax.YTickLabel = arrayfun(@(x) sprintf(axis3_tick_fmt, x), axis3_values, 'UniformOutput', false);
    axis square tight;
    grid on;
    ax.GridAlpha = 0.3;
end

function contour_levels_out = normalize_contour_levels(contour_levels_in)
    contour_levels_out = contour_levels_in(:)';
    if numel(contour_levels_out) == 1
        % MATLAB interprets scalar contour argument as "number of levels".
        % Duplicate the value so it is treated as a specific iso-level.
        contour_levels_out = [contour_levels_out contour_levels_out];
    end
end

function apply_colormap(manual_colormap)
    if ~isempty(manual_colormap)
        if isnumeric(manual_colormap) && size(manual_colormap, 2) == 3
            colormap(manual_colormap);
            return;
        end
        if ischar(manual_colormap) || (isstring(manual_colormap) && isscalar(manual_colormap))
            colormap(char(manual_colormap));
            return;
        end
        error(['manual_colormap must be an N-by-3 numeric colormap matrix ' ...
               'or a valid colormap name string.']);
    end

    try
        colormap(slanCM('viridis'));
    catch
        colormap('parula');
    end
end

function [metric_data, metric_name] = get_metric_data_mode1(dataset_metrics, metric_type, loader_mode)
    if loader_mode == 2
        error('Requested mode=1 plot, but loaded metrics axis mode is 2. Reload with loadMetricDataset_new(1).');
    end
    switch lower(metric_type)
        case 'kernel'
            metric_data = average_over_repetitions(dataset_metrics.kernel_quality_final);
            metric_name = 'Kernel Similarity';
        case 'combined'
            metric_data = average_over_repetitions(dataset_metrics.combined_activationScore);
            metric_name = 'Combined Activation Score';
        case 'fidelity'
            metric_data = average_over_repetitions(get_observation_fidelity(dataset_metrics, 1));
            metric_name = 'Observation Fidelity';
        case 'multiplied'
            k = average_over_repetitions(dataset_metrics.kernel_quality_final);
            a = average_over_repetitions(dataset_metrics.activation_similarity_final);
            metric_data = k .* a;
            metric_name = 'Kernel x Activation';
        case 'kernel_baseline'
            if ~isfield(dataset_metrics, 'kernel_similarity_baseline_final')
                error(['metric_type ''kernel_baseline'' requires ' ...
                       'dataset_metrics.kernel_similarity_baseline_final.']);
            end
            metric_data = average_over_repetitions(dataset_metrics.kernel_similarity_baseline_final);
            metric_name = 'Kernel Similarity Baseline';
        otherwise
            error('metric_type must be ''kernel'', ''combined'', ''fidelity'', ''multiplied'', or ''kernel_baseline''.');
    end
end

function [metric_data, metric_name] = get_metric_data_mode2(dataset_metrics, metric_type, loader_mode)
    if loader_mode == 2
        suffix = '';
    else
        suffix = '_by_side_length_ratio';
    end

    switch lower(metric_type)
        case 'kernel'
            metric_data = average_over_repetitions(select_field(dataset_metrics, ['kernel_quality_final' suffix], 'kernel_quality_final'));
            metric_name = 'Kernel Similarity';
        case 'combined'
            metric_data = average_over_repetitions(select_field(dataset_metrics, ['combined_activationScore' suffix], 'combined_activationScore'));
            metric_name = 'Combined Activation Score';
        case 'fidelity'
            metric_data = average_over_repetitions(get_observation_fidelity(dataset_metrics, 2));
            metric_name = 'Observation Fidelity';
        case 'multiplied'
            k = average_over_repetitions(select_field(dataset_metrics, ['kernel_quality_final' suffix], 'kernel_quality_final'));
            a = average_over_repetitions(select_field(dataset_metrics, ['activation_similarity_final' suffix], 'activation_similarity_final'));
            metric_data = k .* a;
            metric_name = 'Kernel x Activation';
        case 'kernel_baseline'
            metric_data = average_over_repetitions(select_field(dataset_metrics, ...
                ['kernel_similarity_baseline' suffix], 'kernel_similarity_baseline_final'));
            metric_name = 'Kernel Similarity Baseline';
        otherwise
            error('metric_type must be ''kernel'', ''combined'', ''fidelity'', ''multiplied'', or ''kernel_baseline''.');
    end
end

function fidelity_data = get_observation_fidelity(dataset_metrics, requested_mode)
    if isfield(dataset_metrics, 'observation_fidelity_axis_mode')
        stored_mode = dataset_metrics.observation_fidelity_axis_mode;
        if stored_mode ~= requested_mode
            error(['Observation Fidelity is currently stored in axis mode %d, but plot mode %d was requested. ' ...
                   'Rebuild it with build_observation_fidelity_metrics(metrics, %d).'], ...
                   stored_mode, requested_mode, requested_mode);
        end
    end

    if requested_mode == 2
        fidelity_data = select_field(dataset_metrics, 'observation_fidelity_by_side_length_ratio', 'observation_fidelity');
        if isfield(dataset_metrics, 'side_length_ratio_values') && isfield(dataset_metrics, 'Nobs_values')
            n_axis3 = size(fidelity_data, 3);
            n_side = numel(dataset_metrics.side_length_ratio_values);
            n_nobs = numel(dataset_metrics.Nobs_values);
            if n_axis3 == n_nobs && n_axis3 ~= n_side
                error(['Observation Fidelity appears to be stored on Nobs axis (size=%d), but mode=2 was requested. ' ...
                       'Rebuild it with build_observation_fidelity_metrics(metrics, 2).'], n_axis3);
            end
        end
    else
        if ~isfield(dataset_metrics, 'observation_fidelity')
            error('metric_type ''fidelity'' requires dataset_metrics.observation_fidelity.');
        end
        fidelity_data = dataset_metrics.observation_fidelity;
        if isfield(dataset_metrics, 'side_length_ratio_values') && isfield(dataset_metrics, 'Nobs_values')
            n_axis3 = size(fidelity_data, 3);
            n_side = numel(dataset_metrics.side_length_ratio_values);
            n_nobs = numel(dataset_metrics.Nobs_values);
            if n_axis3 == n_side && n_axis3 ~= n_nobs
                error(['Observation Fidelity appears to be stored on side_length_ratio axis (size=%d), but mode=1 was requested. ' ...
                       'Rebuild it with build_observation_fidelity_metrics(metrics, 1).'], n_axis3);
            end
        end
    end
end

function metric_data = align_metric_tensor_layout(metric_data, n_snr, n_theta, n_axis3, metric_type)
    target = [n_snr, n_theta, n_axis3];
    sz = size(metric_data);
    if numel(sz) < 3
        error('Metric data for ''%s'' must be at least 3D ([SNR, theta, axis3]).', metric_type);
    end
    sz3 = sz(1:3);
    if isequal(sz3, target)
        return;
    end

    perms3 = [1 2 3; 1 3 2; 2 1 3; 2 3 1; 3 1 2; 3 2 1];
    matched_perm = [];
    for i = 1:size(perms3, 1)
        if isequal(sz3(perms3(i, :)), target)
            matched_perm = perms3(i, :);
            break;
        end
    end

    if isempty(matched_perm)
        error(['Metric ''%s'' has unexpected size [%d %d %d]; expected a permutation of [%d %d %d] ' ...
               'for [SNR, theta, axis3].'], metric_type, sz3(1), sz3(2), sz3(3), target(1), target(2), target(3));
    end

    if numel(sz) > 3
        tail = 4:numel(sz);
        metric_data = permute(metric_data, [matched_perm, tail]);
    else
        metric_data = permute(metric_data, matched_perm);
    end
end

function avg_data = average_over_repetitions(data)
    if ndims(data) == 4
        avg_data = mean(data, 4, 'omitnan');
    else
        avg_data = data;
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
