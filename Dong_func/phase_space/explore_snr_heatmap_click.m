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
    loader_mode = get_loader_axis_mode(dataset_metrics);

    switch mode
        case 1
            axis3_values = dataset_metrics.Nobs_values(:)';
            axis3_label = 'N_{obs}';
            axis3_tick_fmt = '%.0f';
            [metric_data, metric_name] = get_metric_data_mode1(dataset_metrics, metric_type, loader_mode);
        case 2
            if ~isfield(dataset_metrics, 'side_length_ratio_values') || isempty(dataset_metrics.side_length_ratio_values)
                error('mode=2 requires dataset_metrics.side_length_ratio_values.');
            end
            axis3_values = dataset_metrics.side_length_ratio_values(:)';
            axis3_label = 'side length ratio';
            axis3_tick_fmt = '%.3f';
            [metric_data, metric_name] = get_metric_data_mode2(dataset_metrics, metric_type, loader_mode);
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
            axis_idx = axis3_idx;
        else
            if loader_mode == 2
                axis_idx = axis3_idx;
            else
                axis_idx = resolve_nobs_from_side_ratio(dataset_metrics, snr_idx, theta_idx, axis3_idx);
                if isempty(axis_idx)
                    fprintf('No reconstruction data found for selected side-length-ratio point.\n');
                    continue;
                end
            end
        end

        [rep_idx, rep_val] = resolve_available_rep(dataset_metrics, snr_idx, theta_idx, axis_idx);
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
                click_count, snr_used, theta_pick, read_nobs_for_slot(dataset_metrics, snr_idx, theta_idx, axis_idx, rep_idx), rep_val);
        else
            fprintf('\nClick %d -> nearest grid: SNR=%.2f, density=%.2e, side ratio=%.4f, Nobs=%d, rep=%d\n', ...
                click_count, snr_used, theta_pick, axis3_pick, read_nobs_for_slot(dataset_metrics, snr_idx, theta_idx, axis_idx, rep_idx), rep_val);
        end

        [Y, Y_clean, A0, A0_noiseless, X0, Aout, Xout, bout, extras] = ...
            get_dataset_payload(dataset_metrics, snr_idx, theta_idx, axis_idx, rep_idx);

        if isempty(Aout) || isempty(Xout) || isempty(Y)
            fprintf('Dataset payload is incomplete for this point. Skipping visualization.\n');
            continue;
        end

        if mode == 1
            detail_title = sprintf('#%d | SNR=%.2f, density=%.2e, Nobs=%d, rep=%d', ...
                click_count, snr_used, theta_pick, read_nobs_for_slot(dataset_metrics, snr_idx, theta_idx, axis_idx, rep_idx), rep_val);
        else
            detail_title = sprintf('#%d | SNR=%.2f, density=%.2e, side ratio=%.4f, Nobs=%d, rep=%d', ...
                click_count, snr_used, theta_pick, axis3_pick, read_nobs_for_slot(dataset_metrics, snr_idx, theta_idx, axis_idx, rep_idx), rep_val);
        end
        visualizeResults_click(Y, Aout, X0, Xout, bout, detail_title);
        plot_kernel_comparison_grid(A0, Aout, A0_noiseless, detail_title);
        plot_activation_matching_grid(X0, Xout, Aout, detail_title);

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
        colormap(slanCM('viridis'));
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
    ax.YScale = 'linear';
    ax.YTick = axis3_values;
    apply_defect_density_tick_style(ax, theta_values);
    ax.YTickLabel = arrayfun(@(x) sprintf(axis3_tick_fmt, x), axis3_values, 'UniformOutput', false);
    axis square tight;
    grid on;
    ax.GridAlpha = 0.3;
end

function [metric_data, metric_name] = get_metric_data_mode1(dataset_metrics, metric_type, loader_mode)
    if loader_mode == 2
        error('Requested mode=1 explorer, but loaded metrics axis mode is 2. Reload with loadMetricDataset_new(1).');
    end
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
        case 'multiplied'
            k = average_over_repetitions(select_field(dataset_metrics, ['kernel_quality_final' suffix], 'kernel_quality_final'));
            a = average_over_repetitions(select_field(dataset_metrics, ['activation_similarity_final' suffix], 'activation_similarity_final'));
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

function nobs_val = read_nobs_for_slot(dm, snr_idx, theta_idx, axis_idx, rep_idx)
    if isfield(dm, 'Nobs_at_axis3') && ~isempty(dm.Nobs_at_axis3)
        if ndims(dm.Nobs_at_axis3) == 4
            nobs_val = dm.Nobs_at_axis3(snr_idx, theta_idx, axis_idx, rep_idx);
        else
            nobs_val = dm.Nobs_at_axis3(snr_idx, theta_idx, axis_idx);
        end
        if ~isnan(nobs_val)
            return;
        end
    end
    if isfield(dm, 'Nobs_values') && axis_idx <= numel(dm.Nobs_values)
        nobs_val = dm.Nobs_values(axis_idx);
    else
        nobs_val = NaN;
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
    target_side_ratio = dm.side_length_ratio_values(side_idx);
    for ni = 1:numel(dm.Nobs_values)
        for r = 1:rep_count
            if isempty(dm.Aout{snr_idx, theta_idx, ni, r})
                continue;
            end
            if isfield(dm, 'side_length_ratio_at_axis3') && ~isempty(dm.side_length_ratio_at_axis3)
                ratio_here = dm.side_length_ratio_at_axis3(snr_idx, theta_idx, ni, r);
                if isfinite(ratio_here) && abs(ratio_here - target_side_ratio) < 1e-10
                    nobs_idx = ni;
                    return;
                end
            elseif isfield(dm, 'kernel_quality_final_by_side_length_ratio')
                side_val = dm.kernel_quality_final_by_side_length_ratio(snr_idx, theta_idx, side_idx, r);
                nobs_val = dm.kernel_quality_final(snr_idx, theta_idx, ni, r);
                if ~isnan(side_val) && ~isnan(nobs_val) && abs(side_val - nobs_val) < 1e-10
                    nobs_idx = ni;
                    return;
                end
            else
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
    get_dataset_payload(dm, snr_idx, theta_idx, axis_idx, rep_idx)

    if ndims(dm.Y) < 4
        Y = dm.Y{snr_idx, theta_idx, axis_idx};
        Y_clean = dm.Y_clean{snr_idx, theta_idx, axis_idx};
        A0 = dm.A0{snr_idx, theta_idx, axis_idx};
        A0_noiseless = dm.A0_noiseless{snr_idx, theta_idx, axis_idx};
        X0 = dm.X0{snr_idx, theta_idx, axis_idx};
        Aout = dm.Aout{snr_idx, theta_idx, axis_idx};
        Xout = dm.Xout{snr_idx, theta_idx, axis_idx};
        bout = dm.bout{snr_idx, theta_idx, axis_idx};
        extras = dm.extras{snr_idx, theta_idx, axis_idx};
    else
        Y = dm.Y{snr_idx, theta_idx, axis_idx, rep_idx};
        Y_clean = dm.Y_clean{snr_idx, theta_idx, axis_idx, rep_idx};
        A0 = dm.A0{snr_idx, theta_idx, axis_idx, rep_idx};
        A0_noiseless = dm.A0_noiseless{snr_idx, theta_idx, axis_idx, rep_idx};
        X0 = dm.X0{snr_idx, theta_idx, axis_idx, rep_idx};
        Aout = dm.Aout{snr_idx, theta_idx, axis_idx, rep_idx};
        Xout = dm.Xout{snr_idx, theta_idx, axis_idx, rep_idx};
        bout = dm.bout{snr_idx, theta_idx, axis_idx, rep_idx};
        extras = dm.extras{snr_idx, theta_idx, axis_idx, rep_idx};
    end
end

function loader_mode = get_loader_axis_mode(dm)
    if isfield(dm, 'axis3_mode') && isscalar(dm.axis3_mode)
        loader_mode = dm.axis3_mode;
    else
        loader_mode = 1;
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

function plot_kernel_comparison_grid(A0, Aout, A0_noiseless, detail_title)
    A0_cells = to_kernel_cells(A0);
    Aout_cells = to_kernel_cells(Aout);
    A0n_cells = to_kernel_cells(A0_noiseless);

    n = max([numel(A0_cells), numel(Aout_cells), numel(A0n_cells), 1]);
    fig = figure('Name', sprintf('Kernel panels | %s', detail_title), ...
                 'Position', [140 80 1450 max(320, 260 * n)]);

    col_titles = {'A0', 'Aout', 'A0\_noiseless'};
    for k = 1:n
        kernels_row = {get_kernel_at(A0_cells, k), get_kernel_at(Aout_cells, k), ...
                       get_kernel_at(A0n_cells, k)};
        for c = 1:3
            subplot(n, 3, (k - 1) * 3 + c);
            ker = kernels_row{c};
            if isempty(ker)
                axis off;
                text(0.5, 0.5, 'N/A', 'HorizontalAlignment', 'center', ...
                    'VerticalAlignment', 'middle', 'Color', [0.35 0.35 0.35], 'FontSize', 10);
            else
                imagesc(ker);
                axis image;
                axis off;
                colorbar;
            end
            if k == 1
                title(col_titles{c}, 'Interpreter', 'none');
            end
            if c == 1
                ylabel(sprintf('Kernel %d', k));
            end
        end
    end
    sgtitle(sprintf('Kernel comparison | %s', detail_title), 'Interpreter', 'none');
    set(fig, 'Color', 'w');
end

function cells = to_kernel_cells(v)
    cells = {};
    if isempty(v)
        return;
    end
    if iscell(v)
        if numel(v) == 1 && iscell(v{1})
            cells = v{1};
        else
            cells = v;
        end
        return;
    end
    if isnumeric(v)
        if ndims(v) == 2
            cells = {v};
        elseif ndims(v) >= 3
            cells = cell(1, size(v, 3));
            for i = 1:size(v, 3)
                cells{i} = v(:, :, i);
            end
        end
    end
end

function ker = get_kernel_at(cells, idx)
    ker = [];
    if isempty(cells) || idx > numel(cells)
        return;
    end
    ker = cells{idx};
end

function plot_activation_matching_grid(X0, Xout, Aout, detail_title)
% Show activation-map matching diagnostics using the same pipeline
% used by visualizeResults.m (alignment + adaptive filtering).
    if isempty(X0) || isempty(Xout) || isempty(Aout)
        return;
    end
    if ndims(X0) ~= 3 || ndims(Xout) ~= 3
        return;
    end
    if ~isequal(size(X0), size(Xout))
        fprintf('Activation map plotting skipped: X0/Xout size mismatch.\n');
        return;
    end

    num_kernels = numel(Aout);
    if size(X0, 3) ~= num_kernels
        fprintf('Activation map plotting skipped: kernel/channel mismatch.\n');
        return;
    end

    kernel_size = zeros(num_kernels, 2);
    for k = 1:num_kernels
        ak = Aout{k};
        if isempty(ak) || ~ismatrix(ak)
            fprintf('Activation map plotting skipped: invalid kernel at index %d.\n', k);
            return;
        end
        kernel_size(k, :) = size(ak);
    end

    try
        evaluateActivationReconstruction(X0, Xout, kernel_size, true);
        sgtitle(sprintf('Activation matching | %s', detail_title), 'Interpreter', 'none');
    catch ME
        fprintf('Activation map plotting failed: %s\n', ME.message);
    end
end
