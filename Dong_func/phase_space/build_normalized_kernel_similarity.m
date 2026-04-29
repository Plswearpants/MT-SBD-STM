function metrics = build_normalized_kernel_similarity(metrics, epsilon)
% Build normalized kernel similarity from loaded dataset fields.
%
% This function expects output from loadMetricDataset_new.m and adds:
%   metrics.normalized_kernel_similarity_final
%   metrics.kernel_similarity_baseline_final
% Optional (if side_length_ratio_values exists):
%   metrics.normalized_kernel_similarity_by_side_length_ratio
%   metrics.kernel_similarity_baseline_by_side_length_ratio
%
% Definition per kernel k:
%   s_k = corr2(Q(Aout_k), Q(A0_noiseless_k))
%   b_k = corr2(Q(A0_k), Q(A0_noiseless_k))
%   s_k_norm = s_k / max(b_k, epsilon)
%
% where Q(.) is normalized QPI magnitude:
%   Q(A) = abs(fftshift(fft2(A))) / max(abs(fftshift(fft2(A))), [], 'all')

    if nargin < 2 || isempty(epsilon)
        epsilon = 1e-8;
    end

    required_fields = {'Aout', 'A0', 'A0_noiseless'};
    for i = 1:numel(required_fields)
        if ~isfield(metrics, required_fields{i})
            error('metrics must contain field "%s".', required_fields{i});
        end
    end

    data_dims = size(metrics.Aout);
    num_dims = numel(data_dims);
    if num_dims < 3 || num_dims > 4
        error('metrics.Aout must be indexed as [SNR, density, N_obs, rep] (3D/4D).');
    end

    nks_tensor = nan(data_dims);
    baseline_tensor = nan(data_dims);
    total_slots = numel(metrics.Aout);
    wb = waitbar(0, 'Building normalized kernel similarity...');
    cleanup_obj = onCleanup(@() close_waitbar_safe(wb));

    loader_mode = get_loader_axis_mode(metrics);
    has_side_grid = isfield(metrics, 'side_length_ratio_values') && ~isempty(metrics.side_length_ratio_values) && ...
                    isfield(metrics, 'SNR_values') && isfield(metrics, 'theta_cap_values');
    aggregate_to_side = has_side_grid && loader_mode ~= 2;
    if aggregate_to_side
        side_values = metrics.side_length_ratio_values(:);
        num_snr = numel(metrics.SNR_values);
        num_density = numel(metrics.theta_cap_values);
        if num_dims == 4
            num_rep = data_dims(4);
        else
            num_rep = 1;
        end
        side_sum = zeros(num_snr, num_density, numel(side_values), num_rep);
        side_count = zeros(num_snr, num_density, numel(side_values), num_rep);
        side_baseline_sum = zeros(num_snr, num_density, numel(side_values), num_rep);
        side_baseline_count = zeros(num_snr, num_density, numel(side_values), num_rep);
    end

    for linear_idx = 1:total_slots
        if num_dims == 4
            [idx_snr, idx_density, idx_nobs, idx_rep] = ind2sub(data_dims, linear_idx);
        else
            [idx_snr, idx_density, idx_nobs] = ind2sub(data_dims, linear_idx);
            idx_rep = 1;
        end

        Aout_raw = get_cell_entry(metrics.Aout, num_dims, idx_snr, idx_density, idx_nobs, idx_rep);
        A0_raw = get_cell_entry(metrics.A0, num_dims, idx_snr, idx_density, idx_nobs, idx_rep);
        A0_clean_raw = get_cell_entry(metrics.A0_noiseless, num_dims, idx_snr, idx_density, idx_nobs, idx_rep);

        Aout = to_kernel_cell(Aout_raw);
        A0_noisy = to_kernel_cell(A0_raw);
        A0_clean = to_kernel_cell(A0_clean_raw);

        if isempty(Aout) || isempty(A0_noisy) || isempty(A0_clean)
            update_waitbar(wb, linear_idx, total_slots);
            continue;
        end

        num_kernels = min([numel(Aout), numel(A0_noisy), numel(A0_clean)]);
        if num_kernels < 1
            update_waitbar(wb, linear_idx, total_slots);
            continue;
        end

        per_kernel_nks = nan(1, num_kernels);
        per_kernel_baseline = nan(1, num_kernels);
        for k = 1:num_kernels
            if isempty(Aout{k}) || isempty(A0_noisy{k}) || isempty(A0_clean{k})
                continue;
            end

            s_k = fftspace_corr2(Aout{k}, A0_clean{k});
            b_k = fftspace_corr2(A0_noisy{k}, A0_clean{k});
            denom = max(b_k, epsilon);
            per_kernel_nks(k) = s_k / denom;
            per_kernel_baseline(k) = b_k;
        end

        nks_value = mean(per_kernel_nks, 'omitnan');
        baseline_value = mean(per_kernel_baseline, 'omitnan');
        if isnan(nks_value) && isnan(baseline_value)
            update_waitbar(wb, linear_idx, total_slots);
            continue;
        end

        nks_tensor(linear_idx) = nks_value;
        baseline_tensor(linear_idx) = baseline_value;

        ratio_est = estimate_side_ratio(metrics, num_dims, idx_snr, idx_density, idx_nobs, idx_rep, A0_clean);

        if aggregate_to_side && ~isnan(ratio_est)
            [~, side_idx] = min(abs(side_values - ratio_est));
            if is_ratio_match(ratio_est, side_values(side_idx))
                if ~isnan(nks_value)
                    side_sum(idx_snr, idx_density, side_idx, idx_rep) = ...
                        side_sum(idx_snr, idx_density, side_idx, idx_rep) + nks_value;
                    side_count(idx_snr, idx_density, side_idx, idx_rep) = ...
                        side_count(idx_snr, idx_density, side_idx, idx_rep) + 1;
                end
                if ~isnan(baseline_value)
                    side_baseline_sum(idx_snr, idx_density, side_idx, idx_rep) = ...
                        side_baseline_sum(idx_snr, idx_density, side_idx, idx_rep) + baseline_value;
                    side_baseline_count(idx_snr, idx_density, side_idx, idx_rep) = ...
                        side_baseline_count(idx_snr, idx_density, side_idx, idx_rep) + 1;
                end
            end
        end
        update_waitbar(wb, linear_idx, total_slots);
    end

    metrics.normalized_kernel_similarity_final = nks_tensor;
    metrics.kernel_similarity_baseline_final = baseline_tensor;

    if aggregate_to_side
        nks_side = nan(size(side_sum));
        valid = side_count > 0;
        nks_side(valid) = side_sum(valid) ./ side_count(valid);
        metrics.normalized_kernel_similarity_by_side_length_ratio = nks_side;

        baseline_side = nan(size(side_baseline_sum));
        valid_baseline = side_baseline_count > 0;
        baseline_side(valid_baseline) = side_baseline_sum(valid_baseline) ./ side_baseline_count(valid_baseline);
        metrics.kernel_similarity_baseline_by_side_length_ratio = baseline_side;
    end
    metrics.normalized_kernel_similarity_axis_mode = loader_mode;

    delete(cleanup_obj);
    close_waitbar_safe(wb);
end

function update_waitbar(wb, idx, total_slots)
    if ~ishandle(wb)
        return;
    end
    if idx == 1 || idx == total_slots || mod(idx, max(1, floor(total_slots / 100))) == 0
        waitbar(idx / total_slots, wb, sprintf('Building normalized kernel similarity... %d/%d', idx, total_slots));
    end
end

function close_waitbar_safe(wb)
    if ~isempty(wb) && ishandle(wb)
        close(wb);
    end
end

function kernels = to_kernel_cell(v)
    kernels = {};
    if isempty(v)
        return;
    end

    if iscell(v)
        if numel(v) == 1 && iscell(v{1})
            kernels = v{1};
        else
            kernels = v;
        end
        return;
    end

    if isnumeric(v)
        if ndims(v) == 2
            kernels = {v};
            return;
        end
        if ndims(v) >= 3
            kernels = cell(1, size(v, 3));
            for k = 1:size(v, 3)
                kernels{k} = v(:, :, k);
            end
        end
    end
end

function ratio_est = estimate_side_ratio(metrics, num_dims, idx_snr, idx_density, idx_nobs, idx_rep, clean_kernels)
    ratio_est = NaN;

    kernel_side = get_kernel_side(clean_kernels);
    if isnan(kernel_side)
        return;
    end

    if isfield(metrics, 'Y')
        Y = get_cell_entry(metrics.Y, num_dims, idx_snr, idx_density, idx_nobs, idx_rep);
        if ~isempty(Y)
            y_size = size(Y);
            if numel(y_size) >= 2
                n_obs = max(y_size(1:2));
                if n_obs > 0
                    ratio_est = kernel_side / n_obs;
                    return;
                end
            end
        end
    end

    if isfield(metrics, 'Nobs_values') && idx_nobs <= numel(metrics.Nobs_values)
        n_obs = metrics.Nobs_values(idx_nobs);
        if ~isempty(n_obs) && n_obs > 0
            ratio_est = kernel_side / n_obs;
        end
    end
end

function kernel_side = get_kernel_side(kernels)
    kernel_side = NaN;
    if isempty(kernels)
        return;
    end

    dims = [];
    for k = 1:numel(kernels)
        if isempty(kernels{k})
            continue;
        end
        sz = size(kernels{k});
        if numel(sz) >= 2
            dims(end+1) = max(sz(1:2)); %#ok<AGROW>
        end
    end

    if ~isempty(dims)
        kernel_side = max(dims);
    end
end

function tf = is_ratio_match(ratio_estimate, target_ratio)
    rel_tol = 5e-2;
    tf = abs(ratio_estimate - target_ratio) <= rel_tol * max(abs(target_ratio), eps);
end

function v = get_cell_entry(c, num_dims, idx_snr, idx_density, idx_nobs, idx_rep)
    if num_dims == 3
        v = c{idx_snr, idx_density, idx_nobs};
    else
        v = c{idx_snr, idx_density, idx_nobs, idx_rep};
    end
end

function mode = get_loader_axis_mode(metrics)
    if isfield(metrics, 'axis3_mode') && isscalar(metrics.axis3_mode)
        mode = metrics.axis3_mode;
    else
        mode = 1;
    end
end
