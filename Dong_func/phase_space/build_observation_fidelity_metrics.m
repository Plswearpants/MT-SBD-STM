function metrics = build_observation_fidelity_metrics(metrics, mode)
% Build per-dataset observation fidelity metrics from loaded fields.
%
% Definition per dataset slot:
%   noise    = Y - Y_clean
%   residual = Y - sum_i (convfft2(Aout_i, Xout_i) + b_i)
%   observation_fidelity = std(noise, 0, 'all') / std(residual, 0, 'all')
%
% Adds fields:
%   metrics.observation_fidelity
%   metrics.observation_noise_std
%   metrics.observation_residual_std
%
% Axis mode:
%   mode = 1 (default): save tensors indexed by N_obs
%   mode = 2          : save tensors indexed by side_length_ratio
%
% Notes:
%   - Only one axis version is saved to avoid duplicated storage.
%   - For mode=2, values are aggregated from N_obs slots to side_ratio bins.

    if nargin < 2 || isempty(mode)
        mode = 1;
    end
    if ~(isscalar(mode) && any(mode == [1, 2]))
        error('mode must be 1 (N_obs) or 2 (side_length_ratio).');
    end

    loader_mode = get_loader_axis_mode(metrics);
    if mode == 1 && loader_mode == 2
        error('mode=1 requested, but metrics were loaded in axis mode 2. Reload with loadMetricDataset_new(1).');
    end

    required_fields = {'Y', 'Y_clean', 'Aout', 'Xout'};
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

    fidelity_tensor = nan(data_dims);
    noise_std_tensor = nan(data_dims);
    residual_std_tensor = nan(data_dims);

    total_slots = numel(metrics.Aout);
    wb = waitbar(0, 'Building observation fidelity metrics...');
    cleanup_obj = onCleanup(@() close_waitbar_safe(wb));

    has_side_grid = isfield(metrics, 'side_length_ratio_values') && ~isempty(metrics.side_length_ratio_values) && ...
                    isfield(metrics, 'SNR_values') && isfield(metrics, 'theta_cap_values');
    if mode == 2 && ~has_side_grid
        error('mode=2 requires metrics.side_length_ratio_values, SNR_values, and theta_cap_values.');
    end
    aggregate_to_side = mode == 2 && loader_mode ~= 2;
    if aggregate_to_side
        side_values = metrics.side_length_ratio_values(:);
        num_snr = numel(metrics.SNR_values);
        num_density = numel(metrics.theta_cap_values);
        if num_dims == 4
            num_rep = data_dims(4);
        else
            num_rep = 1;
        end

        sum_fid = zeros(num_snr, num_density, numel(side_values), num_rep);
        sum_noise = zeros(num_snr, num_density, numel(side_values), num_rep);
        sum_res = zeros(num_snr, num_density, numel(side_values), num_rep);

        cnt_fid = zeros(num_snr, num_density, numel(side_values), num_rep);
        cnt_noise = zeros(num_snr, num_density, numel(side_values), num_rep);
        cnt_res = zeros(num_snr, num_density, numel(side_values), num_rep);
        fallback_ratio_count = 0;
    end

    for linear_idx = 1:total_slots
        if num_dims == 4
            [idx_snr, idx_density, idx_nobs, idx_rep] = ind2sub(data_dims, linear_idx);
        else
            [idx_snr, idx_density, idx_nobs] = ind2sub(data_dims, linear_idx);
            idx_rep = 1;
        end

        Y = get_cell_entry(metrics.Y, num_dims, idx_snr, idx_density, idx_nobs, idx_rep);
        Y_clean = get_cell_entry(metrics.Y_clean, num_dims, idx_snr, idx_density, idx_nobs, idx_rep);
        Aout_raw = get_cell_entry(metrics.Aout, num_dims, idx_snr, idx_density, idx_nobs, idx_rep);
        Xout_raw = get_cell_entry(metrics.Xout, num_dims, idx_snr, idx_density, idx_nobs, idx_rep);

        if isempty(Y) || isempty(Y_clean) || isempty(Aout_raw) || isempty(Xout_raw)
            update_waitbar(wb, linear_idx, total_slots);
            continue;
        end

        Aout = to_kernel_cells(Aout_raw);
        Xout = to_activation_tensor(Xout_raw);
        if isempty(Aout) || isempty(Xout)
            update_waitbar(wb, linear_idx, total_slots);
            continue;
        end

        bout = get_bout_for_slot(metrics, num_dims, idx_snr, idx_density, idx_nobs, idx_rep, numel(Aout));
        Y_rec = reconstruct_observation(Aout, Xout, bout, size(Y));
        if isempty(Y_rec)
            update_waitbar(wb, linear_idx, total_slots);
            continue;
        end

        noise = Y - Y_clean;
        residual = Y - Y_rec;

        noise_std = std(noise, 0, 'all');
        residual_std = std(residual, 0, 'all');
        if ~isfinite(noise_std) || ~isfinite(residual_std) || residual_std <= 0
            update_waitbar(wb, linear_idx, total_slots);
            continue;
        end

        fidelity = noise_std / residual_std;
        fidelity_tensor(linear_idx) = fidelity;
        noise_std_tensor(linear_idx) = noise_std;
        residual_std_tensor(linear_idx) = residual_std;

        if aggregate_to_side
            ratio_est = resolve_side_ratio_for_slot(metrics, num_dims, idx_snr, idx_density, idx_nobs, idx_rep, Aout, Y);
            if isnan(ratio_est)
                update_waitbar(wb, linear_idx, total_slots);
                continue;
            end
            if ~has_exact_side_ratio_for_slot(metrics, num_dims, idx_snr, idx_density, idx_nobs, idx_rep)
                fallback_ratio_count = fallback_ratio_count + 1;
            end
            if ~isnan(ratio_est)
                [~, side_idx] = min(abs(side_values - ratio_est));
                % Match loader behavior: always assign to nearest designed side-ratio bin.
                sum_fid(idx_snr, idx_density, side_idx, idx_rep) = sum_fid(idx_snr, idx_density, side_idx, idx_rep) + fidelity;
                sum_noise(idx_snr, idx_density, side_idx, idx_rep) = sum_noise(idx_snr, idx_density, side_idx, idx_rep) + noise_std;
                sum_res(idx_snr, idx_density, side_idx, idx_rep) = sum_res(idx_snr, idx_density, side_idx, idx_rep) + residual_std;
                cnt_fid(idx_snr, idx_density, side_idx, idx_rep) = cnt_fid(idx_snr, idx_density, side_idx, idx_rep) + 1;
                cnt_noise(idx_snr, idx_density, side_idx, idx_rep) = cnt_noise(idx_snr, idx_density, side_idx, idx_rep) + 1;
                cnt_res(idx_snr, idx_density, side_idx, idx_rep) = cnt_res(idx_snr, idx_density, side_idx, idx_rep) + 1;
            end
        end

        update_waitbar(wb, linear_idx, total_slots);
    end

    if mode == 1
        metrics.observation_fidelity = fidelity_tensor;
        metrics.observation_noise_std = noise_std_tensor;
        metrics.observation_residual_std = residual_std_tensor;
    else
        if aggregate_to_side
            metrics.observation_fidelity = divide_with_counts(sum_fid, cnt_fid);
            metrics.observation_noise_std = divide_with_counts(sum_noise, cnt_noise);
            metrics.observation_residual_std = divide_with_counts(sum_res, cnt_res);
            if fallback_ratio_count > 0
                warning('build_observation_fidelity_metrics:SideRatioFallback', ...
                    ['Used geometry-based side-ratio fallback for %d slot(s) because ' ...
                     'exact side_length_ratio_at_axis3 metadata was unavailable.'], ...
                    fallback_ratio_count);
            end
        else
            metrics.observation_fidelity = fidelity_tensor;
            metrics.observation_noise_std = noise_std_tensor;
            metrics.observation_residual_std = residual_std_tensor;
        end
    end

    metrics.observation_fidelity_axis_mode = mode;

    delete(cleanup_obj);
    close_waitbar_safe(wb);
end

function Y_rec = reconstruct_observation(Aout, Xout, bout, y_size)
    num_kernels = min(numel(Aout), size(Xout, 3));
    if num_kernels < 1
        Y_rec = [];
        return;
    end

    Y_rec = zeros(y_size);
    for k = 1:num_kernels
        if isempty(Aout{k})
            continue;
        end
        Y_rec = Y_rec + convfft2(Aout{k}, Xout(:, :, k)) + bout(k);
    end
end

function bout_vec = get_bout_for_slot(metrics, num_dims, idx_snr, idx_density, idx_nobs, idx_rep, num_kernels)
    bout_vec = zeros(num_kernels, 1);
    if ~isfield(metrics, 'bout')
        return;
    end

    bout_raw = get_cell_entry(metrics.bout, num_dims, idx_snr, idx_density, idx_nobs, idx_rep);
    if isempty(bout_raw)
        return;
    end

    bout_raw = bout_raw(:);
    if numel(bout_raw) == 1
        bout_vec = repmat(bout_raw, num_kernels, 1);
    elseif numel(bout_raw) < num_kernels
        bout_vec(1:numel(bout_raw)) = bout_raw;
    else
        bout_vec = bout_raw(1:num_kernels);
    end
end

function arr = divide_with_counts(s, c)
    arr = nan(size(s));
    valid = c > 0;
    arr(valid) = s(valid) ./ c(valid);
end

function update_waitbar(wb, idx, total_slots)
    if ~ishandle(wb)
        return;
    end
    if idx == 1 || idx == total_slots || mod(idx, max(1, floor(total_slots / 100))) == 0
        waitbar(idx / total_slots, wb, sprintf('Building observation fidelity metrics... %d/%d', idx, total_slots));
    end
end

function close_waitbar_safe(wb)
    if ~isempty(wb) && ishandle(wb)
        close(wb);
    end
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

function X = to_activation_tensor(v)
    X = [];
    if isempty(v)
        return;
    end
    if iscell(v)
        if numel(v) == 1 && isnumeric(v{1})
            X = v{1};
        elseif isnumeric(v{1})
            X = cat(3, v{:});
        end
        return;
    end
    if isnumeric(v)
        X = v;
    end
end

function ratio_est = resolve_side_ratio_for_slot(metrics, num_dims, idx_snr, idx_density, idx_nobs, idx_rep, kernels, Y)
    ratio_exact = get_side_ratio_at_axis3_slot(metrics, num_dims, idx_snr, idx_density, idx_nobs, idx_rep);
    if ~isnan(ratio_exact)
        ratio_est = ratio_exact;
        return;
    end

    ratio_est = estimate_side_ratio_from_geometry(metrics, num_dims, idx_nobs, kernels, Y);
end

function tf = has_exact_side_ratio_for_slot(metrics, num_dims, idx_snr, idx_density, idx_nobs, idx_rep)
    tf = ~isnan(get_side_ratio_at_axis3_slot(metrics, num_dims, idx_snr, idx_density, idx_nobs, idx_rep));
end

function ratio = get_side_ratio_at_axis3_slot(metrics, num_dims, idx_snr, idx_density, idx_nobs, idx_rep)
    ratio = NaN;
    if ~isfield(metrics, 'side_length_ratio_at_axis3')
        return;
    end
    try
        if num_dims == 3
            ratio = metrics.side_length_ratio_at_axis3(idx_snr, idx_density, idx_nobs);
        else
            ratio = metrics.side_length_ratio_at_axis3(idx_snr, idx_density, idx_nobs, idx_rep);
        end
    catch
        ratio = NaN;
    end
    if ~isfinite(ratio)
        ratio = NaN;
    end
end

function ratio_est = estimate_side_ratio_from_geometry(metrics, num_dims, idx_nobs, kernels, Y)
    ratio_est = NaN;
    kernel_side = get_kernel_side(kernels);
    if isnan(kernel_side)
        return;
    end

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
            dims(end + 1) = max(sz(1:2)); %#ok<AGROW>
        end
    end
    if ~isempty(dims)
        kernel_side = max(dims);
    end
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
