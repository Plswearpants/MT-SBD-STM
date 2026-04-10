function explore_snr_heatmap_click(dataset_metrics, snr_value, metric_type, mode, interp_factor)
% Interactive explorer on a fixed-SNR interpolated heatmap.
% User clicks points, nearest dataset is selected and visualized.
%
% Inputs:
%   dataset_metrics : struct from loadMetricDataset_new
%   snr_value       : target SNR value (nearest grid value is used)
%   metric_type     : 'kernel' | 'combined' | 'multiplied' (default 'kernel')
%   mode            : 1 -> axis-3 is Nobs, 2 -> axis-3 is side_length_ratio (default 1)
%   interp_factor   : interpolation upsampling factor (default 5)

    if nargin < 2 || isempty(snr_value)
        error('Please provide snr_value.');
    end
    if nargin < 3 || isempty(metric_type)
        metric_type = 'kernel';
    end
    if nargin < 4 || isempty(mode)
        mode = 1;
    end
    if nargin < 5 || isempty(interp_factor)
        interp_factor = 5;
    end

    switch mode
        case 1
            axis3_values = dataset_metrics.Nobs_values(:)';
            axis3_label = 'N_{obs}';
            axis3_tick_fmt = '%.0f';
            [metric_data, metric_name] = get_metric_data_mode1(dataset_metrics, metric_type);
        case 2
            if ~isfield(dataset_metrics, 'side_length_ratio_values') || isempty(dataset_metrics.side_length_ratio_values)
                error('mode=2 requires dataset_metrics.side_length_ratio_values.');
            end
            axis3_values = dataset_metrics.side_length_ratio_values(:)';
            axis3_label = 'side length ratio';
            axis3_tick_fmt = '%.3f';
            [metric_data, metric_name] = get_metric_data_mode2(dataset_metrics, metric_type);
        otherwise
            error('mode must be 1 (Nobs) or 2 (side_length_ratio).');
    end

    theta_values = dataset_metrics.theta_cap_values(:)';
    snr_values = dataset_metrics.SNR_values(:)';
    [~, snr_idx] = min(abs(snr_values - snr_value));
    snr_used = snr_values(snr_idx);

    metric_slice = squeeze(metric_data(snr_idx, :, :)); % [theta x axis3]

    h_heatmap_fig = figure('Name', sprintf('Interactive %s at SNR=%.2f', metric_name, snr_used), ...
           'Position', [120 120 760 640]);
    plot_interpolated_slice(theta_values, axis3_values, metric_slice, ...
        sprintf('%s | SNR=%.2f', metric_name, snr_used), axis3_label, axis3_tick_fmt, interp_factor);
    h_heatmap_ax = gca;
    hold on;

    fprintf('\nInteractive selection mode:\n');
    fprintf('- Left click: choose nearest dataset point on parameter grid.\n');
    fprintf('- After each selection, confirm whether to continue clicking.\n');
    fprintf('- Press Enter or right click at any click prompt to finish.\n');

    click_colors = lines(20);
    click_count = 0;
    while true
        [x_click, y_click, button] = ginput(1);
        if isempty(button) || button ~= 1
            break;
        end
        if x_click <= 0
            fprintf('Ignored click with non-positive x on log axis.\n');
            continue;
        end

        click_count = click_count + 1;
        color_now = click_colors(mod(click_count - 1, size(click_colors, 1)) + 1, :);

        [~, theta_idx] = min(abs(log10(theta_values) - log10(x_click)));
        [~, axis3_idx] = min(abs(axis3_values - y_click));

        if mode == 1
            nobs_idx = axis3_idx;
        else
            nobs_idx = resolve_nobs_from_side_ratio(dataset_metrics, snr_idx, theta_idx, axis3_idx);
            if isempty(nobs_idx)
                fprintf('No reconstruction data found for selected side-length-ratio point.\n');
                continue;
            end
        end

        [rep_idx, rep_val] = resolve_available_rep(dataset_metrics, snr_idx, theta_idx, nobs_idx);
        if isempty(rep_idx)
            fprintf('No reconstruction data found at nearest grid point.\n');
            continue;
        end

        theta_pick = theta_values(theta_idx);
        axis3_pick = axis3_values(axis3_idx);
        plot(theta_pick, axis3_pick, 'o', 'Color', color_now, 'MarkerSize', 11, 'LineWidth', 2);
        text(theta_pick, axis3_pick, sprintf(' #%d', click_count), ...
            'Color', color_now, 'FontWeight', 'bold', 'FontSize', 10, ...
            'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');

        if mode == 1
            fprintf('\nClick %d -> nearest grid: SNR=%.2f, density=%.2e, Nobs=%d, rep=%d\n', ...
                click_count, snr_used, theta_pick, dataset_metrics.Nobs_values(nobs_idx), rep_val);
        else
            fprintf('\nClick %d -> nearest grid: SNR=%.2f, density=%.2e, side ratio=%.4f, Nobs=%d, rep=%d\n', ...
                click_count, snr_used, theta_pick, axis3_pick, dataset_metrics.Nobs_values(nobs_idx), rep_val);
        end

        [Y, Y_clean, A0, A0_noiseless, X0, Aout, Xout, bout, extras] = ...
            get_dataset_payload(dataset_metrics, snr_idx, theta_idx, nobs_idx, rep_idx);

        if isempty(Aout) || isempty(Xout) || isempty(Y)
            fprintf('Dataset payload is incomplete for this point. Skipping visualization.\n');
            continue;
        end

        if mode == 1
            detail_title = sprintf('#%d | SNR=%.2f, density=%.2e, Nobs=%d, rep=%d', ...
                click_count, snr_used, theta_pick, dataset_metrics.Nobs_values(nobs_idx), rep_val);
        else
            detail_title = sprintf('#%d | SNR=%.2f, density=%.2e, side ratio=%.4f, Nobs=%d, rep=%d', ...
                click_count, snr_used, theta_pick, axis3_pick, dataset_metrics.Nobs_values(nobs_idx), rep_val);
        end
        visualizeResults_click(Y, Aout, X0, Xout, bout, detail_title);

        continue_choice = input('Do another click? (y/n): ', 's');
        if isempty(continue_choice) || ~strcmpi(strtrim(continue_choice), 'y')
            break;
        end
        if isgraphics(h_heatmap_fig)
            figure(h_heatmap_fig);
            if isgraphics(h_heatmap_ax)
                axes(h_heatmap_ax);
            end
            drawnow;
        end
    end

    hold off;
end

function plot_interpolated_slice(theta_values, axis3_values, metric_slice, title_str, axis3_label, axis3_tick_fmt, interp_factor)
    theta_log = log10(theta_values(:)');
    y_vals = axis3_values(:)';

    [X, Y] = meshgrid(theta_log, y_vals);
    Z = metric_slice';

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

    xlabel('defect density');
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

function [rep_idx, rep_val] = resolve_available_rep(dm, snr_idx, theta_idx, nobs_idx)
    rep_idx = [];
    rep_val = [];

    if ndims(dm.Aout) < 4
        if ~isempty(dm.Aout{snr_idx, theta_idx, nobs_idx})
            rep_idx = 1;
            rep_val = 1;
        end
        return;
    end

    rep_values = 1:size(dm.Aout, 4);
    if isfield(dm, 'repetition_values') && numel(dm.repetition_values) == numel(rep_values)
        rep_values = dm.repetition_values;
    end

    for r = 1:numel(rep_values)
        if ~isempty(dm.Aout{snr_idx, theta_idx, nobs_idx, r}) && ...
           ~isempty(dm.Xout{snr_idx, theta_idx, nobs_idx, r})
            rep_idx = r;
            rep_val = rep_values(r);
            return;
        end
    end
end

function nobs_idx = resolve_nobs_from_side_ratio(dm, snr_idx, theta_idx, side_idx)
% Find an N_obs index with reconstruction payload that matches this side ratio.
    nobs_idx = [];

    if ndims(dm.Aout) < 4
        for ni = 1:numel(dm.Nobs_values)
            if ~isempty(dm.Aout{snr_idx, theta_idx, ni})
                nobs_idx = ni;
                return;
            end
        end
        return;
    end

    rep_count = size(dm.Aout, 4);
    for ni = 1:numel(dm.Nobs_values)
        for r = 1:rep_count
            if isempty(dm.Aout{snr_idx, theta_idx, ni, r})
                continue;
            end
            side_val = dm.kernel_quality_final_by_side_length_ratio(snr_idx, theta_idx, side_idx, r);
            nobs_val = dm.kernel_quality_final(snr_idx, theta_idx, ni, r);
            if ~isnan(side_val) && ~isnan(nobs_val) && abs(side_val - nobs_val) < 1e-10
                nobs_idx = ni;
                return;
            end
        end
    end

    % Fallback: first available N_obs at this (SNR, theta)
    for ni = 1:numel(dm.Nobs_values)
        for r = 1:rep_count
            if ~isempty(dm.Aout{snr_idx, theta_idx, ni, r})
                nobs_idx = ni;
                return;
            end
        end
    end
end

function [Y, Y_clean, A0, A0_noiseless, X0, Aout, Xout, bout, extras] = ...
    get_dataset_payload(dm, snr_idx, theta_idx, nobs_idx, rep_idx)

    if ndims(dm.Y) < 4
        Y = dm.Y{snr_idx, theta_idx, nobs_idx};
        Y_clean = dm.Y_clean{snr_idx, theta_idx, nobs_idx};
        A0 = dm.A0{snr_idx, theta_idx, nobs_idx};
        A0_noiseless = dm.A0_noiseless{snr_idx, theta_idx, nobs_idx};
        X0 = dm.X0{snr_idx, theta_idx, nobs_idx};
        Aout = dm.Aout{snr_idx, theta_idx, nobs_idx};
        Xout = dm.Xout{snr_idx, theta_idx, nobs_idx};
        bout = dm.bout{snr_idx, theta_idx, nobs_idx};
        extras = dm.extras{snr_idx, theta_idx, nobs_idx};
    else
        Y = dm.Y{snr_idx, theta_idx, nobs_idx, rep_idx};
        Y_clean = dm.Y_clean{snr_idx, theta_idx, nobs_idx, rep_idx};
        A0 = dm.A0{snr_idx, theta_idx, nobs_idx, rep_idx};
        A0_noiseless = dm.A0_noiseless{snr_idx, theta_idx, nobs_idx, rep_idx};
        X0 = dm.X0{snr_idx, theta_idx, nobs_idx, rep_idx};
        Aout = dm.Aout{snr_idx, theta_idx, nobs_idx, rep_idx};
        Xout = dm.Xout{snr_idx, theta_idx, nobs_idx, rep_idx};
        bout = dm.bout{snr_idx, theta_idx, nobs_idx, rep_idx};
        extras = dm.extras{snr_idx, theta_idx, nobs_idx, rep_idx};
    end
end
