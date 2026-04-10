function metrics = build_nor_loo_metrics(metrics)
% Build metric-style NOR/LOO tensors from existing loaded dataset fields.
%
% Adds:
%   metrics.NOR
%   metrics.LOO
%   metrics.side_length_ratio_estimate
% Optional (if side_length_ratio_values exists):
%   metrics.NOR_by_side_length_ratio
%   metrics.LOO_by_side_length_ratio

    if ~isfield(metrics, 'X0') || ~isfield(metrics, 'Y')
        error('metrics must contain X0 and Y fields.');
    end
    if ~isfield(metrics, 'SNR_values') || ~isfield(metrics, 'theta_cap_values')
        error('metrics must contain SNR_values and theta_cap_values.');
    end

    x0_dims = size(metrics.X0);
    num_dims = numel(x0_dims);
    if num_dims < 3 || num_dims > 4
        error('metrics.X0 must be 3D or 4D indexed as [SNR, density, N_obs, rep].');
    end

    nor_tensor = nan(x0_dims);
    loo_tensor = nan(x0_dims);
    side_ratio_tensor = nan(x0_dims);

    total_slots = numel(metrics.X0);
    wb = waitbar(0, 'Building NOR/LOO metrics...');
    cleanup_obj = onCleanup(@() close_waitbar_safe(wb));

    for linear_idx = 1:total_slots
        if num_dims == 3
            [idx_snr, idx_density, idx_nobs] = ind2sub(x0_dims, linear_idx);
            idx_rep = 1;
        else
            [idx_snr, idx_density, idx_nobs, idx_rep] = ind2sub(x0_dims, linear_idx);
        end

        X0 = get_cell_entry(metrics.X0, num_dims, idx_snr, idx_density, idx_nobs, idx_rep);
        if isempty(X0)
            update_waitbar(wb, linear_idx, total_slots);
            continue;
        end

        Y = get_cell_entry(metrics.Y, num_dims, idx_snr, idx_density, idx_nobs, idx_rep);
        if isempty(Y)
            y_size = size(X0(:, :, 1));
        else
            y_size = size(Y);
            y_size = y_size(1:2);
        end

        kernel_sizes = infer_kernel_sizes(metrics, num_dims, idx_snr, idx_density, idx_nobs, idx_rep, size(X0, 3));
        if isempty(kernel_sizes)
            update_waitbar(wb, linear_idx, total_slots);
            continue;
        end

        [nor_val, loo_val] = compute_nor_loo_from_x0(X0, y_size, kernel_sizes);
        nor_tensor(linear_idx) = nor_val;
        loo_tensor(linear_idx) = loo_val;
        side_ratio_tensor(linear_idx) = estimate_side_length_ratio(y_size, kernel_sizes);

        update_waitbar(wb, linear_idx, total_slots);
    end

    % Keep computed tensors in metric-style indexing
    metrics.NOR = nor_tensor;
    metrics.LOO = loo_tensor;
    metrics.side_length_ratio_estimate = side_ratio_tensor;

    if isfield(metrics, 'side_length_ratio_values') && ~isempty(metrics.side_length_ratio_values)
        [nor_side, loo_side] = aggregate_by_side_ratio(metrics, nor_tensor, loo_tensor, side_ratio_tensor, num_dims, x0_dims);
        metrics.NOR_by_side_length_ratio = nor_side;
        metrics.LOO_by_side_length_ratio = loo_side;
    end

    delete(cleanup_obj);
    close_waitbar_safe(wb);
end

function update_waitbar(wb, idx, total_slots)
    if ~ishandle(wb)
        return;
    end
    if idx == 1 || idx == total_slots || mod(idx, max(1, floor(total_slots / 100))) == 0
        waitbar(idx / total_slots, wb, sprintf('Building NOR/LOO metrics... %d/%d', idx, total_slots));
    end
end

function close_waitbar_safe(wb)
    if ~isempty(wb) && ishandle(wb)
        close(wb);
    end
end

function kernel_sizes = infer_kernel_sizes(metrics, num_dims, idx_snr, idx_density, idx_nobs, idx_rep, num_kernels)
    kernel_sizes = [];

    if isfield(metrics, 'A0')
        A0 = get_cell_entry(metrics.A0, num_dims, idx_snr, idx_density, idx_nobs, idx_rep);
        kernel_sizes = kernel_sizes_from_A0(A0, num_kernels);
        if ~isempty(kernel_sizes)
            return;
        end
    end

    if isfield(metrics, 'A0_noiseless')
        A0n = get_cell_entry(metrics.A0_noiseless, num_dims, idx_snr, idx_density, idx_nobs, idx_rep);
        kernel_sizes = kernel_sizes_from_A0(A0n, num_kernels);
    end
end

function kernel_sizes = kernel_sizes_from_A0(A0, num_kernels)
    kernel_sizes = [];
    if isempty(A0)
        return;
    end

    if iscell(A0)
        n = min(num_kernels, numel(A0));
        kernel_sizes = zeros(num_kernels, 2);
        for k = 1:n
            if isempty(A0{k})
                continue;
            end
            sz = size(A0{k});
            kernel_sizes(k, :) = sz(1:2);
        end
        if any(kernel_sizes(:, 1) == 0)
            valid = kernel_sizes(:, 1) > 0;
            if any(valid)
                fallback = kernel_sizes(find(valid, 1), :); %#ok<FNDSB>
                kernel_sizes(~valid, :) = repmat(fallback, sum(~valid), 1);
            else
                kernel_sizes = [];
            end
        end
    elseif isnumeric(A0) && ndims(A0) >= 2
        sz = size(A0);
        if ndims(A0) == 2
            kernel_sizes = repmat(sz(1:2), num_kernels, 1);
        else
            n = min(num_kernels, sz(3));
            kernel_sizes = zeros(num_kernels, 2);
            for k = 1:n
                kernel_sizes(k, :) = sz(1:2);
            end
            if n < num_kernels && n >= 1
                kernel_sizes(n+1:end, :) = repmat(kernel_sizes(1, :), num_kernels - n, 1);
            end
        end
    end
end

function [nor_val, loo_val] = compute_nor_loo_from_x0(X0, y_size, kernel_sizes)
    H = y_size(1);
    W = y_size(2);
    coverage_map = zeros(H, W);

    num_kernels = size(X0, 3);
    for k = 1:num_kernels
        act_map = X0(:, :, k) > 0;
        [rows, cols] = find(act_map);
        if isempty(rows)
            continue;
        end

        ksize = kernel_sizes(min(k, size(kernel_sizes, 1)), :);
        for j = 1:numel(rows)
            [r0, r1, c0, c1] = kernel_patch_bounds(rows(j), cols(j), ksize, H, W);
            coverage_map(r0:r1, c0:c1) = coverage_map(r0:r1, c0:c1) + 1;
        end
    end

    total_area = H * W;
    nor_val = nnz(coverage_map == 1) / total_area;
    loo_val = sum(coverage_map, 'all') / total_area;
end

function [r0, r1, c0, c1] = kernel_patch_bounds(center_row, center_col, ksize, H, W)
    half_h = floor(ksize(1) / 2);
    half_w = floor(ksize(2) / 2);

    r_desired = (center_row - half_h):(center_row + half_h);
    c_desired = (center_col - half_w):(center_col + half_w);

    if numel(r_desired) > ksize(1)
        r_desired = r_desired(1:ksize(1));
    end
    if numel(c_desired) > ksize(2)
        c_desired = c_desired(1:ksize(2));
    end

    r0 = max(1, min(r_desired));
    r1 = min(H, max(r_desired));
    c0 = max(1, min(c_desired));
    c1 = min(W, max(c_desired));
end

function ratio_est = estimate_side_length_ratio(y_size, kernel_sizes)
    M_max = max(kernel_sizes(:));
    N_obs_est = max(y_size(1:2));
    ratio_est = M_max / N_obs_est;
end

function [nor_side, loo_side] = aggregate_by_side_ratio(metrics, nor_tensor, loo_tensor, ratio_tensor, num_dims, x0_dims)
    side_values = metrics.side_length_ratio_values(:);
    num_snr = numel(metrics.SNR_values);
    num_density = numel(metrics.theta_cap_values);
    num_side = numel(side_values);
    if num_dims == 4
        num_rep = x0_dims(4);
    else
        num_rep = 1;
    end

    nor_side = nan(num_snr, num_density, num_side, num_rep);
    loo_side = nan(num_snr, num_density, num_side, num_rep);

    for s = 1:num_snr
        for d = 1:num_density
            for r = 1:num_rep
                if num_dims == 4
                    nor_slice = squeeze(nor_tensor(s, d, :, r));
                    loo_slice = squeeze(loo_tensor(s, d, :, r));
                    ratio_slice = squeeze(ratio_tensor(s, d, :, r));
                else
                    nor_slice = squeeze(nor_tensor(s, d, :));
                    loo_slice = squeeze(loo_tensor(s, d, :));
                    ratio_slice = squeeze(ratio_tensor(s, d, :));
                end

                for si = 1:num_side
                    target_ratio = side_values(si);
                    ratio_mask = is_ratio_match(ratio_slice, target_ratio);
                    if ~any(ratio_mask)
                        continue;
                    end
                    nor_side(s, d, si, r) = mean(nor_slice(ratio_mask), 'omitnan');
                    loo_side(s, d, si, r) = mean(loo_slice(ratio_mask), 'omitnan');
                end
            end
        end
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
