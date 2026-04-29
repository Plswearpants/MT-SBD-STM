function [side_ratio_points, occurrence_points, kernel_similarity_points, combined_activation_points] = plot_occurrence_vs_length_ratio_by_snr(metrics, snr_value)
%PLOT_OCCURRENCE_VS_LENGTH_RATIO_BY_SNR
% Plot average defect occurrence vs side-length ratio for one SNR.
%
% Inputs:
%   metrics   : dataset_metrics returned by loadMetricDataset_new
%   snr_value : optional numeric SNR value; if omitted, prompts user
%
% Outputs:
%   side_ratio_points : side-length ratio value per scatter point
%   occurrence_points : average occurrence value per scatter point
%   kernel_similarity_points   : kernel quality value per point
%   combined_activation_points : combined activation value per point

    if nargin < 2
        snr_value = [];
    end

    if ~isfield(metrics, 'SNR_values') || isempty(metrics.SNR_values)
        error('metrics.SNR_values is required.');
    end
    if ~isfield(metrics, 'X0') || isempty(metrics.X0)
        error('metrics.X0 is required to compute occurrence counts.');
    end
    if ~isfield(metrics, 'kernel_quality_final') || isempty(metrics.kernel_quality_final)
        error('metrics.kernel_quality_final is required for colormap.');
    end
    if ~isfield(metrics, 'combined_activationScore') || isempty(metrics.combined_activationScore)
        error('metrics.combined_activationScore is required for colormap.');
    end

    if isempty(snr_value)
        fprintf('\nAvailable SNR range: [%.3g, %.3g]\n', ...
            min(metrics.SNR_values), max(metrics.SNR_values));
        snr_value = input('Enter SNR value to plot: ');
    end

    [~, snr_idx] = min(abs(metrics.SNR_values - snr_value));
    snr_actual = metrics.SNR_values(snr_idx);
    fprintf('Using nearest SNR on grid: %.6g\n', snr_actual);

    dims_x0 = size(metrics.X0);
    if numel(dims_x0) < 3
        error('metrics.X0 must be indexed as [SNR, theta, N_obs] or [SNR, theta, N_obs, rep].');
    end
    num_theta = dims_x0(2);
    num_nobs = dims_x0(3);
    if numel(dims_x0) >= 4
        num_rep = dims_x0(4);
    else
        num_rep = 1;
    end
    has_rep = ndims(metrics.X0) == 4;

    max_points = num_theta * num_nobs * num_rep;
    side_ratio_points = nan(max_points, 1);
    occurrence_points = nan(max_points, 1);
    kernel_similarity_points = nan(max_points, 1);
    combined_activation_points = nan(max_points, 1);
    point_count = 0;

    for theta_idx = 1:num_theta
        for nobs_idx = 1:num_nobs
            for rep_idx = 1:num_rep
                if has_rep
                    X0 = metrics.X0{snr_idx, theta_idx, nobs_idx, rep_idx};
                    Y = get_cell_if_present(metrics, 'Y', snr_idx, theta_idx, nobs_idx, rep_idx);
                    A0 = get_cell_if_present(metrics, 'A0_noiseless', snr_idx, theta_idx, nobs_idx, rep_idx);
                else
                    X0 = metrics.X0{snr_idx, theta_idx, nobs_idx};
                    Y = get_cell_if_present(metrics, 'Y', snr_idx, theta_idx, nobs_idx);
                    A0 = get_cell_if_present(metrics, 'A0_noiseless', snr_idx, theta_idx, nobs_idx);
                end

                if isempty(X0)
                    continue;
                end

                num_kernels = size(X0, 3);
                if num_kernels == 0
                    continue;
                end
                occ_counts = zeros(1, num_kernels);
                for k = 1:num_kernels
                    occ_counts(k) = sum(X0(:,:,k), 'all');
                end
                avg_occurrence = mean(occ_counts);

                ratio_val = estimate_side_length_ratio_from_entry(metrics, Y, A0, nobs_idx);
                if isnan(ratio_val)
                    continue;
                end

                if has_rep
                    kernel_val = metrics.kernel_quality_final(snr_idx, theta_idx, nobs_idx, rep_idx);
                    combined_val = metrics.combined_activationScore(snr_idx, theta_idx, nobs_idx, rep_idx);
                else
                    kernel_val = metrics.kernel_quality_final(snr_idx, theta_idx, nobs_idx);
                    combined_val = metrics.combined_activationScore(snr_idx, theta_idx, nobs_idx);
                end

                point_count = point_count + 1;
                side_ratio_points(point_count) = ratio_val;
                occurrence_points(point_count) = avg_occurrence;
                kernel_similarity_points(point_count) = kernel_val;
                combined_activation_points(point_count) = combined_val;
            end
        end
    end

    side_ratio_points = side_ratio_points(1:point_count);
    occurrence_points = occurrence_points(1:point_count);
    kernel_similarity_points = kernel_similarity_points(1:point_count);
    combined_activation_points = combined_activation_points(1:point_count);

    if isempty(side_ratio_points)
        warning('No valid points found for SNR=%.6g.', snr_actual);
        return;
    end

    % Log axis requires strictly positive occurrence values.
    positive_mask = occurrence_points > 0;
    side_ratio_points = side_ratio_points(positive_mask);
    occurrence_points = occurrence_points(positive_mask);
    kernel_similarity_points = kernel_similarity_points(positive_mask);
    combined_activation_points = combined_activation_points(positive_mask);
    if isempty(occurrence_points)
        warning('All occurrence points are non-positive for SNR=%.6g; cannot use log Y axis.', snr_actual);
        return;
    end

    try
        custom_cmap = slanCM('viridis');
    catch
        custom_cmap = parula(256);
    end

    figure('Name', sprintf('Occurrence vs side-length ratio (SNR=%.6g)', snr_actual), ...
           'Position', [100 100 1300 500]);
    t = tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

    nexttile;
    mask_kernel = ~isnan(kernel_similarity_points);
    scatter(side_ratio_points(mask_kernel), occurrence_points(mask_kernel), 55, ...
        kernel_similarity_points(mask_kernel), 'filled', ...
        'MarkerFaceAlpha', 0.85, 'MarkerEdgeAlpha', 0.85);
    xlabel('Side-length ratio');
    ylabel('Average occurrence of defects');
    set(gca, 'YScale', 'log');
    title('Colored by kernel similarity');
    grid on;
    colormap(gca, custom_cmap);
    cb1 = colorbar;
    cb1.Label.String = 'Kernel similarity';

    nexttile;
    mask_combined = ~isnan(combined_activation_points);
    scatter(side_ratio_points(mask_combined), occurrence_points(mask_combined), 55, ...
        combined_activation_points(mask_combined), 'filled', ...
        'MarkerFaceAlpha', 0.85, 'MarkerEdgeAlpha', 0.85);
    xlabel('Side-length ratio');
    ylabel('Average occurrence of defects');
    set(gca, 'YScale', 'log');
    title('Colored by combined activation score');
    grid on;
    colormap(gca, custom_cmap);
    cb2 = colorbar;
    cb2.Label.String = 'Combined activation score';

    title(t, sprintf('Occurrence vs side-length ratio (SNR = %.6g)', snr_actual));
end

function value = get_cell_if_present(metrics, field_name, i1, i2, i3, i4)
    value = [];
    if ~isfield(metrics, field_name)
        return;
    end
    arr = metrics.(field_name);
    dims_arr = size(arr);
    if nargin < 7
        if numel(dims_arr) < 3 || any([i1, i2, i3] > dims_arr(1:3))
            return;
        end
        value = arr{i1, i2, i3};
    else
        if numel(dims_arr) >= 4
            if any([i1, i2, i3, i4] > dims_arr(1:4))
                return;
            end
            value = arr{i1, i2, i3, i4};
        elseif any([i1, i2, i3] > dims_arr(1:3))
            return;
        else
            value = arr{i1, i2, i3};
        end
    end
end

function ratio_val = estimate_side_length_ratio_from_entry(metrics, Y, A0, nobs_idx)
    ratio_val = nan;

    if ~isempty(Y) && ~isempty(A0)
        max_kernel_side = infer_max_kernel_side(A0);
        if ~isnan(max_kernel_side)
            n_obs_est = max(size(Y, 1), size(Y, 2));
            if n_obs_est > 0
                ratio_raw = max_kernel_side / n_obs_est;
                ratio_val = snap_to_ratio_grid_if_available(metrics, ratio_raw);
                return;
            end
        end
    end

    if isfield(metrics, 'Nobs_values') && nobs_idx <= numel(metrics.Nobs_values) ...
            && ~isempty(A0)
        max_kernel_side = infer_max_kernel_side(A0);
        if ~isnan(max_kernel_side) && metrics.Nobs_values(nobs_idx) > 0
            ratio_raw = max_kernel_side / metrics.Nobs_values(nobs_idx);
            ratio_val = snap_to_ratio_grid_if_available(metrics, ratio_raw);
        end
    end
end

function max_side = infer_max_kernel_side(A0)
    max_side = nan;
    if isempty(A0)
        return;
    end
    if iscell(A0)
        kernel_sides = nan(numel(A0), 1);
        for i = 1:numel(A0)
            if ~isempty(A0{i})
                kernel_sides(i) = max(size(A0{i}));
            end
        end
        if any(~isnan(kernel_sides))
            max_side = max(kernel_sides, [], 'omitnan');
        end
    else
        max_side = max(size(A0));
    end
end

function snapped = snap_to_ratio_grid_if_available(metrics, ratio_raw)
    snapped = ratio_raw;
    if isfield(metrics, 'side_length_ratio_values') && ~isempty(metrics.side_length_ratio_values)
        [~, idx] = min(abs(metrics.side_length_ratio_values - ratio_raw));
        snapped = metrics.side_length_ratio_values(idx);
    end
end
