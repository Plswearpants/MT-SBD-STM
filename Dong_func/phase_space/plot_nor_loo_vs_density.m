function [density_values, nor_values, loo_values, nor_loo_product] = plot_nor_loo_vs_density(metrics, snr_value, side_length_ratio_value)
% Pure plotting for NOR/LOO and NOR*LOO curves vs defect density.
% Requires precomputed metrics.NOR and metrics.LOO (use build_nor_loo_metrics first).

    if nargin < 2
        snr_value = [];
    end
    if nargin < 3
        side_length_ratio_value = [];
    end

    if ~isfield(metrics, 'NOR') || ~isfield(metrics, 'LOO')
        error('metrics.NOR and metrics.LOO are missing. Run build_nor_loo_metrics(metrics) first.');
    end
    if ~isfield(metrics, 'theta_cap_values')
        error('metrics.theta_cap_values is required for density axis.');
    end
    if ~isfield(metrics, 'SNR_values')
        error('metrics.SNR_values is required for SNR filtering.');
    end

    nor_tensor = metrics.NOR;
    loo_tensor = metrics.LOO;
    x_dims = size(nor_tensor);
    num_dims = numel(x_dims);
    if num_dims < 3 || num_dims > 4
        error('metrics.NOR must be 3D or 4D indexed as [SNR, density, N_obs, rep].');
    end
    if ~isequal(size(nor_tensor), size(loo_tensor))
        error('metrics.NOR and metrics.LOO must have the same size.');
    end
    if any(isnan(nor_tensor(:)))
        warning('plot_nor_loo_vs_density:NORHasNaN', ...
            'metrics.NOR contains NaN values (%d/%d). Some density points may be missing.', ...
            nnz(isnan(nor_tensor)), numel(nor_tensor));
    end
    if any(isnan(loo_tensor(:)))
        warning('plot_nor_loo_vs_density:LOOHasNaN', ...
            'metrics.LOO contains NaN values (%d/%d). Some density points may be missing.', ...
            nnz(isnan(loo_tensor)), numel(loo_tensor));
    end

    if isfield(metrics, 'side_length_ratio_estimate')
        side_ratio_tensor = metrics.side_length_ratio_estimate;
    else
        side_ratio_tensor = nan(x_dims);
    end

    density_values = metrics.theta_cap_values(:);
    num_density = numel(density_values);
    snr_indices = resolve_snr_indices(metrics, snr_value);

    nor_values = nan(num_density, 1);
    loo_values = nan(num_density, 1);
    nor_loo_product = nan(num_density, 1);
    selected_ratio_by_density = nan(num_density, 1);
    for d = 1:num_density
        [nor_vals_d, loo_vals_d, ratio_vals_d, selected_ratio_d] = collect_density_slice( ...
            nor_tensor, loo_tensor, side_ratio_tensor, x_dims, num_dims, snr_indices, d, side_length_ratio_value);
        if ~isempty(nor_vals_d)
            nor_values(d) = mean(nor_vals_d, 'omitnan');
        end
        if ~isempty(loo_vals_d)
            loo_values(d) = mean(loo_vals_d, 'omitnan');
        end
        if ~isempty(nor_vals_d) && ~isempty(loo_vals_d)
            nor_loo_product(d) = mean(nor_vals_d .* loo_vals_d, 'omitnan');
        end
        selected_ratio_by_density(d) = selected_ratio_d;

        % Print actual side-length ratio values used at this density point.
        if ~isempty(ratio_vals_d)
            fprintf('[Density %.6e] actual side-length ratios used: ', density_values(d));
            fprintf('%.6f ', ratio_vals_d);
            fprintf('\n');
        else
            fprintf('[Density %.6e] no valid data points after filtering.\n', density_values(d));
        end
    end

    if ~isempty(side_length_ratio_value)
        valid_sel = selected_ratio_by_density(~isnan(selected_ratio_by_density));
        if ~isempty(valid_sel)
            fprintf('Requested side-length ratio %.6f, closest selected ratios across densities: ', side_length_ratio_value);
            fprintf('%.6f ', unique(valid_sel));
            fprintf('\n');
        end
    end

    if any(isnan(nor_values))
        warning('plot_nor_loo_vs_density:NORCurveHasNaN', ...
            'NOR curve contains NaN at %d/%d density points after filtering.', ...
            nnz(isnan(nor_values)), numel(nor_values));
    end
    if any(isnan(loo_values))
        warning('plot_nor_loo_vs_density:LOOCurveHasNaN', ...
            'LOO curve contains NaN at %d/%d density points after filtering.', ...
            nnz(isnan(loo_values)), numel(loo_values));
    end
    if any(isnan(nor_loo_product))
        warning('plot_nor_loo_vs_density:ProductCurveHasNaN', ...
            'NOR*LOO curve contains NaN at %d/%d density points after filtering.', ...
            nnz(isnan(nor_loo_product)), numel(nor_loo_product));
    end

    figure('Name', 'NOR, LOO, and NOR*LOO vs Density', 'Position', [100 100 1500 450]);
    t = tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

    nexttile;
    plot(density_values, nor_values, '-o', 'LineWidth', 1.8, 'MarkerSize', 6);
    set(gca, 'XScale', 'log');
    xlabel('Defect density');
    ylabel('NOR');
    title('Non-overlapping Ratio vs Density');
    grid on;

    nexttile;
    plot(density_values, loo_values, '-s', 'LineWidth', 1.8, 'MarkerSize', 6);
    set(gca, 'XScale', 'log');
    xlabel('Defect density');
    ylabel('LOO');
    title('Level of Overlapping vs Density');
    grid on;

    nexttile;
    plot(density_values, nor_loo_product, '-d', 'LineWidth', 1.8, 'MarkerSize', 6);
    set(gca, 'XScale', 'log');
    xlabel('Defect density');
    ylabel('NOR * LOO');
    title('NOR*LOO vs Density');
    grid on;

    title(t, build_filter_subtitle(metrics, snr_value, side_length_ratio_value));
end

function snr_indices = resolve_snr_indices(metrics, snr_value)
    if isempty(snr_value)
        snr_indices = 1:numel(metrics.SNR_values);
        return;
    end
    snr_indices = find(abs(metrics.SNR_values - snr_value) < 1e-10, 1);
    if isempty(snr_indices)
        error('Requested SNR %.6g not found in metrics.SNR_values.', snr_value);
    end
end

function [nor_vals_d, loo_vals_d, ratio_vals_d, selected_ratio_d] = collect_density_slice(nor_tensor, loo_tensor, ratio_tensor, x_dims, num_dims, snr_indices, density_idx, side_length_ratio_value)
    nor_vals_d = [];
    loo_vals_d = [];
    ratio_vals_d = [];
    selected_ratio_d = nan;

    total_slots = numel(nor_tensor);
    for linear_idx = 1:total_slots
        if num_dims == 3
            [idx_snr, idx_density, ~] = ind2sub(x_dims, linear_idx);
        else
            [idx_snr, idx_density, ~, ~] = ind2sub(x_dims, linear_idx);
        end

        if ~ismember(idx_snr, snr_indices) || idx_density ~= density_idx
            continue;
        end

        nor_val = nor_tensor(linear_idx);
        loo_val = loo_tensor(linear_idx);
        ratio_est = ratio_tensor(linear_idx);
        if isnan(nor_val) || isnan(loo_val) || isnan(ratio_est)
            continue;
        end

        nor_vals_d(end+1, 1) = nor_val; %#ok<AGROW>
        loo_vals_d(end+1, 1) = loo_val; %#ok<AGROW>
        ratio_vals_d(end+1, 1) = ratio_est; %#ok<AGROW>
    end

    % Closest-match side-length-ratio filtering (replaces tolerance gating).
    if ~isempty(side_length_ratio_value) && ~isempty(ratio_vals_d)
        diffs = abs(ratio_vals_d - side_length_ratio_value);
        min_diff = min(diffs);
        closest_mask = diffs <= (min_diff + eps(max(1, side_length_ratio_value)));
        nor_vals_d = nor_vals_d(closest_mask);
        loo_vals_d = loo_vals_d(closest_mask);
        ratio_vals_d = ratio_vals_d(closest_mask);
        if ~isempty(ratio_vals_d)
            selected_ratio_d = mean(ratio_vals_d, 'omitnan');
        end
    elseif isempty(side_length_ratio_value) && ~isempty(ratio_vals_d)
        selected_ratio_d = mean(ratio_vals_d, 'omitnan');
    end
end

function txt = build_filter_subtitle(metrics, snr_value, side_length_ratio_value)
    if isempty(snr_value)
        snr_txt = 'all SNR';
    else
        snr_txt = sprintf('SNR=%.3g', snr_value);
    end

    if isempty(side_length_ratio_value)
        side_txt = 'all side-length ratios';
    else
        side_txt = sprintf('side-length ratio=%.3g', side_length_ratio_value);
    end

    txt = sprintf('NOR/LOO vs defect density (%s, %s, %d density levels)', ...
        snr_txt, side_txt, numel(metrics.theta_cap_values));
end
