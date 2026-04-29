function h = metrics2heat_general(dataset_metrics, mode, varargin)
% Plot a general 3D heatspace for a user-selected metric.
%
% Inputs:
%   dataset_metrics : struct from loadMetricDataset_new (or compatible)
%   mode            : 1 -> axis-3 is N_obs
%                     2 -> axis-3 is side_length_ratio
%   varargin        : optional name-value settings:
%       'plot_mode'       -> '3d' (default) or 'line_profile'
%       'metric_type'     -> for line_profile: 'combined' or 'kernel'
%       'snr_value'       -> requested SNR value for line_profile
%       'interp_factor'   -> interpolation factor for line_profile (default 4)
%       'manual_colormap' -> colormap matrix or name string
%
% Behavior:
%   - 3D mode: prompts user to select which metric to plot.
%   - Averages over repetition dimension if data is 4D.
%   - line_profile mode: shows one interpolated phase map at chosen SNR and
%     lets user draw a line to inspect metric values along that path.
%
% Output:
%   h : 3D mode -> scatter3 handle (empty if canceled)
%       line_profile mode -> struct with plot handles/profile data

    if nargin < 2 || isempty(mode)
        mode = 1;
    end
    mode = normalize_mode(mode);
    opts = parse_plot_options(varargin{:});
    h = [];

    loader_mode = get_loader_axis_mode(dataset_metrics);
    [axis3_values, axis3_label, axis3_tick_fmt, suffix] = get_axis3_config(dataset_metrics, mode, loader_mode);
    if strcmp(opts.plot_mode, '3d')
        [metric_specs, metric_names] = get_available_metric_options(dataset_metrics, suffix);

        if isempty(metric_specs)
            error('No plottable metric found for mode %d.', mode);
        end

        selected_idx = ask_user_metric_choice(metric_names);
        if isempty(selected_idx)
            fprintf('Metric selection canceled. No plot generated.\n');
            return;
        end

        selected_spec = metric_specs(selected_idx);
        metric_data = evaluate_metric(dataset_metrics, selected_spec, suffix, mode);
        metric_avg = average_over_repetitions(metric_data);

        title_str = sprintf('%s (Averaged over Repetitions)', selected_spec.display_name);
        h = plot_3D_heatspace_general(dataset_metrics.SNR_values, ...
                                      dataset_metrics.theta_cap_values, ...
                                      axis3_values, ...
                                      metric_avg, ...
                                      title_str, ...
                                      axis3_label, ...
                                      axis3_tick_fmt);
    else
        [metric_data, metric_display_name] = get_line_profile_metric(dataset_metrics, suffix, opts.metric_type);
        metric_avg = average_over_repetitions(metric_data);
        [snr_idx, snr_value] = resolve_snr_choice(dataset_metrics.SNR_values, opts.snr_value);
        metric_slice = squeeze(metric_avg(snr_idx, :, :)); % [theta x axis3]

        title_str = sprintf('%s | SNR=%.2f', metric_display_name, snr_value);
        h = plot_phase_space_line_profile(dataset_metrics.theta_cap_values, ...
                                          axis3_values, ...
                                          metric_slice, ...
                                          title_str, ...
                                          axis3_label, ...
                                          axis3_tick_fmt, ...
                                          opts.interp_factor, ...
                                          opts.manual_colormap);
    end
end

function opts = parse_plot_options(varargin)
    p = inputParser;
    p.FunctionName = 'metrics2heat_general';
    addParameter(p, 'plot_mode', '3d', @(x) ischar(x) || (isstring(x) && isscalar(x)));
    addParameter(p, 'metric_type', '', @(x) ischar(x) || (isstring(x) && isscalar(x)));
    addParameter(p, 'snr_value', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x)));
    addParameter(p, 'interp_factor', 4, @(x) isnumeric(x) && isscalar(x) && x >= 1);
    addParameter(p, 'manual_colormap', [], @(x) isempty(x) || isnumeric(x) || ischar(x) || isstring(x));
    parse(p, varargin{:});

    opts = p.Results;
    opts.plot_mode = lower(string(opts.plot_mode));
    if ~any(opts.plot_mode == ["3d", "line_profile"])
        error('plot_mode must be ''3d'' or ''line_profile''.');
    end
    if ~isempty(opts.metric_type)
        opts.metric_type = lower(string(opts.metric_type));
    end
end

function mode_num = normalize_mode(mode)
    if isstring(mode) || ischar(mode)
        mode_num = str2double(string(mode));
    else
        mode_num = mode;
    end

    if ~(isscalar(mode_num) && any(mode_num == [1, 2]))
        error('mode must be 1 (N_obs) or 2 (side_length_ratio).');
    end
end

function [axis3_values, axis3_label, axis3_tick_fmt, suffix] = get_axis3_config(dataset_metrics, mode, loader_mode)
    switch mode
        case 1
            if loader_mode == 2
                error('Requested mode=1 plot, but loaded metrics axis mode is 2. Reload with loadMetricDataset_new(1).');
            end
            if ~isfield(dataset_metrics, 'Nobs_values') || isempty(dataset_metrics.Nobs_values)
                error('mode=1 requires dataset_metrics.Nobs_values.');
            end
            axis3_values = dataset_metrics.Nobs_values;
            axis3_label = 'N_{obs}';
            axis3_tick_fmt = '%.0f';
            suffix = '';
        case 2
            if ~isfield(dataset_metrics, 'side_length_ratio_values') || isempty(dataset_metrics.side_length_ratio_values)
                error('mode=2 requires dataset_metrics.side_length_ratio_values.');
            end
            axis3_values = dataset_metrics.side_length_ratio_values;
            axis3_label = 'side length ratio';
            axis3_tick_fmt = '%.3f';
            if loader_mode == 2
                suffix = '';
            else
                suffix = '_by_side_length_ratio';
            end
    end
end

function [metric_specs, metric_names] = get_available_metric_options(dataset_metrics, suffix)
    metric_specs = struct('id', {}, 'display_name', {}, 'kind', {}, 'field', {});

    metric_specs = append_field_metric(metric_specs, dataset_metrics, ...
        'kernel', 'Kernel Similarity', ['kernel_quality_final' suffix]);
    metric_specs = append_field_metric(metric_specs, dataset_metrics, ...
        'activation', 'Activation Similarity', ['activation_similarity_final' suffix]);
    metric_specs = append_field_metric(metric_specs, dataset_metrics, ...
        'combined', 'Combined Activation Score', ['combined_activationScore' suffix]);
    metric_specs = append_field_metric(metric_specs, dataset_metrics, ...
        'demixing', 'Demixing Score', ['demixing_score' suffix]);
    metric_specs = append_field_metric(metric_specs, dataset_metrics, ...
        'observation_fidelity', 'Observation Fidelity', 'observation_fidelity');
    metric_specs = append_first_existing_metric(metric_specs, dataset_metrics, ...
        'normalized_kernel', 'Normalized Kernel Similarity', ...
        {['normalized_kernel_similarity_final' suffix], ['normalized_kernel_similarity' suffix]});
    metric_specs = append_first_existing_metric(metric_specs, dataset_metrics, ...
        'kernel_baseline', 'Kernel Similarity Baseline', ...
        {['kernel_similarity_baseline_final' suffix], ['kernel_similarity_baseline' suffix]});

    kernel_field = ['kernel_quality_final' suffix];
    activation_field = ['activation_similarity_final' suffix];
    if isfield(dataset_metrics, kernel_field) && isfield(dataset_metrics, activation_field)
        metric_specs(end + 1) = struct( ... %#ok<AGROW>
            'id', 'multiplied', ...
            'display_name', 'Kernel x Activation', ...
            'kind', 'computed', ...
            'field', '');
    end

    metric_names = {metric_specs.display_name};
end

function metric_specs = append_field_metric(metric_specs, dataset_metrics, id, display_name, field_name)
    if isfield(dataset_metrics, field_name)
        metric_specs(end + 1) = struct( ... %#ok<AGROW>
            'id', id, ...
            'display_name', display_name, ...
            'kind', 'field', ...
            'field', field_name);
    end
end

function metric_specs = append_first_existing_metric(metric_specs, dataset_metrics, id, display_name, field_candidates)
    for i = 1:numel(field_candidates)
        field_name = field_candidates{i};
        if isfield(dataset_metrics, field_name)
            metric_specs(end + 1) = struct( ... %#ok<AGROW>
                'id', id, ...
                'display_name', display_name, ...
                'kind', 'field', ...
                'field', field_name);
            return;
        end
    end
end

function selected_idx = ask_user_metric_choice(metric_names)
    selected_idx = [];

    if usejava('desktop')
        [idx, ok] = listdlg( ...
            'PromptString', 'Select metric to plot:', ...
            'SelectionMode', 'single', ...
            'ListString', metric_names, ...
            'ListSize', [320 180]);
        if ok && ~isempty(idx)
            selected_idx = idx;
        end
        return;
    end

    fprintf('\nSelect metric to plot:\n');
    for i = 1:numel(metric_names)
        fprintf('  %d) %s\n', i, metric_names{i});
    end
    idx = input('Enter metric number (empty to cancel): ');
    if isempty(idx)
        return;
    end
    if ~(isscalar(idx) && idx >= 1 && idx <= numel(metric_names) && floor(idx) == idx)
        error('Invalid metric selection index.');
    end
    selected_idx = idx;
end

function metric_data = evaluate_metric(dataset_metrics, metric_spec, suffix, requested_mode)
    switch metric_spec.kind
        case 'field'
            if strcmp(metric_spec.id, 'observation_fidelity') && isfield(dataset_metrics, 'observation_fidelity_axis_mode')
                stored_mode = dataset_metrics.observation_fidelity_axis_mode;
                if stored_mode ~= requested_mode
                    error(['Observation Fidelity is currently stored in axis mode %d, but plot mode %d was requested. ' ...
                           'Rebuild it with build_observation_fidelity_metrics(metrics, %d).'], ...
                           stored_mode, requested_mode, requested_mode);
                end
            end
            metric_data = dataset_metrics.(metric_spec.field);
        case 'computed'
            switch metric_spec.id
                case 'multiplied'
                    metric_data = dataset_metrics.(['kernel_quality_final' suffix]) .* ...
                                  dataset_metrics.(['activation_similarity_final' suffix]);
                otherwise
                    error('Unsupported computed metric id: %s', metric_spec.id);
            end
        otherwise
            error('Unknown metric kind: %s', metric_spec.kind);
    end
end

function [metric_data, metric_display_name] = get_line_profile_metric(dataset_metrics, suffix, requested_metric_type)
    available = struct('id', {}, 'display_name', {}, 'field_name', {});
    available = append_line_metric_if_exists(available, dataset_metrics, ...
        'combined', 'Combined Activation Score', ...
        ['combined_activationScore' suffix], 'combined_activationScore');
    available = append_line_metric_if_exists(available, dataset_metrics, ...
        'kernel', 'Kernel Similarity', ...
        ['kernel_quality_final' suffix], 'kernel_quality_final');

    if isempty(available)
        error('No ''combined'' or ''kernel'' metric found for line_profile mode.');
    end

    selected_idx = [];
    if ~isempty(requested_metric_type)
        for i = 1:numel(available)
            if string(available(i).id) == requested_metric_type
                selected_idx = i;
                break;
            end
        end
        if isempty(selected_idx)
            error('metric_type must be ''combined'' or ''kernel'' (if provided).');
        end
    else
        metric_names = {available.display_name};
        selected_idx = ask_user_metric_choice(metric_names);
        if isempty(selected_idx)
            error('Line-profile metric selection canceled.');
        end
    end

    metric_display_name = available(selected_idx).display_name;
    metric_data = dataset_metrics.(available(selected_idx).field_name);
end

function available = append_line_metric_if_exists(available, dataset_metrics, id, display_name, preferred_field, fallback_field)
    if isfield(dataset_metrics, preferred_field)
        field_name = preferred_field;
    elseif isfield(dataset_metrics, fallback_field)
        field_name = fallback_field;
    else
        return;
    end

    available(end + 1) = struct( ... %#ok<AGROW>
        'id', id, ...
        'display_name', display_name, ...
        'field_name', field_name);
end

function [snr_idx, snr_value] = resolve_snr_choice(snr_values, requested_snr)
    if ~isempty(requested_snr)
        [~, snr_idx] = min(abs(snr_values - requested_snr));
        snr_value = snr_values(snr_idx);
        return;
    end

    snr_labels = arrayfun(@(v) sprintf('SNR = %.2f', v), snr_values, 'UniformOutput', false);
    selected_idx = ask_user_metric_choice(snr_labels);
    if isempty(selected_idx)
        error('SNR selection canceled.');
    end
    snr_idx = selected_idx;
    snr_value = snr_values(snr_idx);
end

function out = plot_phase_space_line_profile(theta_values, axis3_values, metric_slice, title_str, axis3_label, axis3_tick_fmt, interp_factor, manual_colormap)
    theta_log = log10(theta_values(:)');
    y_vals = axis3_values(:)';
    [X, Y] = meshgrid(theta_log, y_vals);
    Z = metric_slice'; % [axis3 x theta]

    xq = linspace(min(theta_log), max(theta_log), max(2, interp_factor * numel(theta_values)));
    yq = linspace(min(y_vals), max(y_vals), max(2, interp_factor * numel(y_vals)));
    [Xq, Yq] = meshgrid(xq, yq);
    Zq = interp2(X, Y, Z, Xq, Yq, 'linear');
    x_plot = 10.^xq;

    fig = figure('Position', [100 100 1200 480]);
    ax_map = subplot(1, 2, 1);
    img_h = imagesc(x_plot, yq, Zq);
    set(ax_map, 'YDir', 'normal');
    axis(ax_map, 'square');
    axis(ax_map, 'tight');
    apply_colormap(manual_colormap);
    clim(ax_map, [0 1]);
    cb = colorbar(ax_map);
    ylabel(cb, 'Performance Score');
    xlabel(ax_map, 'defect density (per site)');
    ylabel(ax_map, axis3_label);
    title(ax_map, sprintf('%s (draw line for profile)', title_str), 'Interpreter', 'none');
    grid(ax_map, 'on');
    ax_map.GridAlpha = 0.3;
    ax_map.YTick = axis3_values;
    apply_defect_density_tick_style(ax_map, theta_values);
    ax_map.YTickLabel = arrayfun(@(v) sprintf(axis3_tick_fmt, v), axis3_values, 'UniformOutput', false);

    ax_profile = subplot(1, 2, 2);
    profile_h = plot(ax_profile, NaN, NaN, 'k-', 'LineWidth', 1.6);
    xlabel(ax_profile, 'Distance along line (log-theta / axis-3 units)');
    ylabel(ax_profile, 'Metric value');
    title(ax_profile, 'Line Profile (live)');
    grid(ax_profile, 'on');

    line_h = gobjects(0);
    live_listeners = [];
    if usejava('desktop')
        try
            x_mid = x_plot(round((numel(x_plot) + 1) / 2));
            y_mid = yq(round((numel(yq) + 1) / 2));
            dx_idx = max(1, floor(numel(x_plot) / 8));
            dy = (max(yq) - min(yq)) * 0.2;
            x1 = x_plot(max(1, round((numel(x_plot) + 1) / 2) - dx_idx));
            x2 = x_plot(min(numel(x_plot), round((numel(x_plot) + 1) / 2) + dx_idx));
            y1 = y_mid - dy;
            y2 = y_mid + dy;
            y1 = min(max(y1, min(yq)), max(yq));
            y2 = min(max(y2, min(yq)), max(yq));

            line_h = drawline(ax_map, ...
                'Color', 'w', ...
                'LineWidth', 1.6, ...
                'Position', [x1, y1; x2, y2]);

            update_profile_plot_from_line(line_h, x_plot, yq, Zq, profile_h, ax_profile);
            live_listeners(1) = addlistener(line_h, 'MovingROI', ...
                @(src, ~) update_profile_plot_from_line(src, x_plot, yq, Zq, profile_h, ax_profile));
            live_listeners(2) = addlistener(line_h, 'ROIMoved', ...
                @(src, ~) update_profile_plot_from_line(src, x_plot, yq, Zq, profile_h, ax_profile));
            wait(line_h);
        catch
            % Fall back to manual two-point selection if drawline is unavailable.
        end
    end

    if isempty(line_h) || ~isvalid(line_h)
        [x_line, y_line, line_h] = request_line_on_axis(ax_map);
        [xs, ys, distance, profile_vals] = compute_line_profile_samples(x_line, y_line, x_plot, yq, Zq);
        set(profile_h, 'XData', distance, 'YData', profile_vals);
        if all(~isnan(profile_vals)) && min(profile_vals) >= 0 && max(profile_vals) <= 1
            ylim(ax_profile, [0 1]);
        end
    else
        line_pos = line_h.Position;
        x_line = line_pos([1, end], 1);
        y_line = line_pos([1, end], 2);
        [xs, ys, distance, profile_vals] = compute_line_profile_samples(x_line, y_line, x_plot, yq, Zq);
    end

    out = struct();
    out.figure = fig;
    out.map_axis = ax_map;
    out.profile_axis = ax_profile;
    out.image_handle = img_h;
    out.drawn_line_handle = line_h;
    out.profile_handle = profile_h;
    out.x_samples = xs;
    out.y_samples = ys;
    out.distance = distance;
    out.profile_values = profile_vals;
    out.live_listeners = live_listeners;
end

function [x_line, y_line, line_h] = request_line_on_axis(ax)
    line_h = gobjects(0);

    if usejava('desktop')
        try
            line_h = drawline(ax, 'Color', 'w', 'LineWidth', 1.6);
            wait(line_h);
            pts = line_h.Position;
            if size(pts, 1) >= 2
                x_line = pts([1, end], 1);
                y_line = pts([1, end], 2);
                return;
            end
        catch
            % Fall back to ginput if drawline is unavailable.
        end
    end

    fprintf('Select two points on the phase map for line profile...\n');
    [x_pick, y_pick] = ginput(2);
    if numel(x_pick) < 2 || numel(y_pick) < 2
        error('Line selection canceled.');
    end
    x_line = x_pick(:)';
    y_line = y_pick(:)';
    hold(ax, 'on');
    line_h = plot(ax, x_line, y_line, 'w--', 'LineWidth', 1.6);
    hold(ax, 'off');
end

function update_profile_plot_from_line(line_roi, x_plot, yq, Zq, profile_h, ax_profile)
    line_pos = line_roi.Position;
    if size(line_pos, 1) < 2
        return;
    end
    x_line = line_pos([1, end], 1);
    y_line = line_pos([1, end], 2);
    [~, ~, distance, profile_vals] = compute_line_profile_samples(x_line, y_line, x_plot, yq, Zq);
    set(profile_h, 'XData', distance, 'YData', profile_vals);
    if all(~isnan(profile_vals)) && min(profile_vals) >= 0 && max(profile_vals) <= 1
        ylim(ax_profile, [0 1]);
    end
    drawnow limitrate;
end

function [xs, ys, distance, profile_vals] = compute_line_profile_samples(x_line, y_line, x_plot, yq, Zq)
    line_log_span = abs(log10(max(x_line) / min(x_line)));
    yq_span = max(yq) - min(yq);
    line_y_span_norm = abs(y_line(2) - y_line(1)) / max(eps, yq_span);
    n_samples = max(300, ceil(1200 * hypot(line_log_span, line_y_span_norm)));

    xs = linspace(x_line(1), x_line(2), n_samples);
    ys = linspace(y_line(1), y_line(2), n_samples);
    profile_vals = interp2(x_plot, yq, Zq, xs, ys, 'linear', NaN);
    distance = [0, cumsum(hypot(diff(log10(xs)), diff(ys)))];
end

function h = plot_3D_heatspace_general(snr_values, theta_values, axis3_values, metric_values, title_str, axis3_label, axis3_tick_fmt)
    [snr_grid, theta_grid, axis3_grid] = meshgrid(snr_values, theta_values, axis3_values);
    snr_grid = permute(snr_grid, [2, 1, 3]);
    theta_grid = permute(theta_grid, [2, 1, 3]);
    axis3_grid = permute(axis3_grid, [2, 1, 3]);

    z_pos = zeros(size(snr_grid));
    for i = 1:numel(snr_values)
        z_pos(snr_grid == snr_values(i)) = i;
    end

    x = theta_grid(:);
    y = axis3_grid(:);
    z = z_pos(:);
    c = metric_values(:);

    valid = ~isnan(c);
    x = x(valid);
    y = y(valid);
    z = z(valid);
    c = c(valid);

    figure('Position', [100 100 900 700]);
    h = scatter3(x, y, z, 90, c, 'filled');

    try
        colormap(slanCM('viridis'));
    catch
        colormap('parula');
    end

    cb = colorbar;
    ylabel(cb, 'Performance Score');

    if ~isempty(c)
        cmin = min(c);
        cmax = max(c);
        if cmin >= 0 && cmax <= 1
            clim([0 1]);
        elseif cmin < cmax
            clim([cmin cmax]);
        end
    end

    xlabel('defect density');
    ylabel(axis3_label);
    zlabel('SNR');
    title(title_str);

    grid on;
    view(-35, 85);

    ax = gca;
    ax.YScale = 'linear';
    ax.ZScale = 'linear';

    ax.YTick = axis3_values;
    ax.ZTick = 1:numel(snr_values);

    apply_defect_density_tick_style(ax, theta_values);
    ax.YTickLabel = arrayfun(@(v) sprintf(axis3_tick_fmt, v), axis3_values, 'UniformOutput', false);
    ax.ZTickLabel = arrayfun(@(v) sprintf('%.2f', v), snr_values, 'UniformOutput', false);

    set(ax, 'XDir', 'normal');
    set(ax, 'YDir', 'reverse');
    set(ax, 'ZDir', 'reverse');
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

function avg_data = average_over_repetitions(data)
    if ndims(data) == 4
        avg_data = mean(data, 4, 'omitnan');
    else
        avg_data = data;
    end
end

function loader_mode = get_loader_axis_mode(dataset_metrics)
    if isfield(dataset_metrics, 'axis3_mode') && isscalar(dataset_metrics.axis3_mode)
        loader_mode = dataset_metrics.axis3_mode;
    else
        loader_mode = 1;
    end
end
