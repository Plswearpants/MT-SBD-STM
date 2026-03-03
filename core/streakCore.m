function varargout = streakCore(op, varargin)
%STREAKCORE Shared logic for streak detection and correction (single source of truth).
%
%   Use one of:
%   [L_for_mask, slider_range] = streakCore('laplacian', data, mode, max_streak_width)
%   mask = streakCore('mask', L_for_mask, mode, threshold)
%   corrected = streakCore('correct', data, mask, method)
%   [var_list, lower_list] = streakCore('sweep', data, mode, max_streak_width, n_thresholds)
%   lower_list = streakCore('threshold_list', L_for_mask, n_thresholds)
%
%   threshold_list: thresholds at linearly spaced *population* ranks (sort L, take 1st, (M/n)th, ...).
%   So resolution follows the data distribution instead of linear value spacing.
%   mode: 'valley' | 'plateau' | 'both'.
%   method: 'neighbor_interp' (column expected from neighbors) | 'background_mean' (replace with non-streak mean) | 'horizontal_avg' (pixel = avg of left/right neighbors).
%
%   'horizontal_avg' only corrects pixels that have both horizontal neighbors (valid_streaks).

    switch lower(op)
        case 'laplacian'
            [varargout{1}, varargout{2}] = computeLaplacian(varargin{1}, varargin{2}, varargin{3});
        case 'mask'
            varargout{1} = maskFromThreshold(varargin{1}, varargin{2}, varargin{3});
        case 'correct'
            varargout{1} = applyCorrection(varargin{1}, varargin{2}, varargin{3});
        case 'sweep'
            [varargout{1}, varargout{2}] = varianceSweep(varargin{1}, varargin{2}, varargin{3}, varargin{4});
        case 'threshold_list'
            varargout{1} = thresholdListFromPopulation(varargin{1}, varargin{2});
        otherwise
            error('streakCore: unknown op "%s".', op);
    end
end

%% -------------------------------------------------------------------------
function [L_for_mask, slider_range] = computeLaplacian(data, mode, max_streak_width)
% Multi-scale 1D Laplacian (x-direction). Same formula everywhere: L_n = left + right - 2*data.
    [rows, cols] = size(data);
    L_combined = zeros(size(data));
    L_mag_combined = zeros(size(data));

    for n = 1:max_streak_width
        data_left_n  = [zeros(rows, n), data(:, 1:end-n)];
        data_right_n = [data(:, n+1:end), zeros(rows, n)];
        L_n = data_left_n + data_right_n - 2*data;
        if strcmpi(mode, 'valley') || strcmpi(mode, 'plateau')
            L_combined = max(L_combined, L_n);
        else
            L_mag_combined = max(L_mag_combined, abs(L_n));
        end
    end

    if strcmpi(mode, 'valley')
        L_for_mask = L_combined;
        mn = 0;
        mx = max(L_for_mask(:));
    elseif strcmpi(mode, 'plateau')
        L_for_mask = L_combined;
        mn = min(L_for_mask(:));
        mx = 0;
        if mn > mx
            [mn, mx] = deal(mx, mn);
        end
    else
        L_for_mask = L_mag_combined;
        mn = min(L_for_mask(:));
        mx = prctile(L_for_mask(:), 99.5);
    end
    if mn == mx
        mx = mn + eps;
    end
    slider_range = [mn, mx];
end

%% -------------------------------------------------------------------------
function mask = maskFromThreshold(L_for_mask, mode, threshold)
    if strcmpi(mode, 'valley')
        mask = (L_for_mask > threshold);
    elseif strcmpi(mode, 'plateau')
        mask = (L_for_mask < threshold);
    else
        mask = (L_for_mask >= threshold);
    end
end

%% -------------------------------------------------------------------------
function corrected = applyCorrection(data, streak_mask, method)
% method: 'neighbor_interp' | 'background_mean' | 'horizontal_avg'
    corrected = data;
    [streak_rows, streak_cols] = find(streak_mask);
    global_mean = mean(data(:));

    if strcmpi(method, 'horizontal_avg')
        % Only pixels that have both left and right neighbors (and exclude col 1 and end)
        valid = streak_mask;
        valid(:,1) = false;
        valid(:,end) = false;
        valid = valid & circshift(streak_mask, [0 -1]) & circshift(streak_mask, [0 1]);
        corrected(valid) = (data(circshift(valid, [0 -1])) + data(circshift(valid, [0 1]))) / 2;
        return;
    end

    unique_cols = unique(streak_cols);
    if strcmpi(method, 'background_mean')
        if any(~streak_mask(:))
            ref_mean = mean(data(~streak_mask));
        else
            ref_mean = global_mean;
        end
    end

    for col = unique_cols'
        rows_in_col = streak_rows(streak_cols == col);
        if isempty(rows_in_col)
            continue;
        end
        streak_avg = mean(data(rows_in_col, col));
        if strcmpi(method, 'background_mean')
            expected_avg = ref_mean;
        else
            % neighbor_interp: expected from left/right column means
            left_col = col - 1;
            right_col = col + 1;
            left_vals = [];
            right_vals = [];
            if left_col >= 1
                left_vals = mean(data(rows_in_col, left_col));
            end
            if right_col <= size(data, 2)
                right_vals = mean(data(rows_in_col, right_col));
            end
            if ~isempty(left_vals) && ~isempty(right_vals)
                expected_avg = left_vals + (col - left_col) * (right_vals - left_vals) / (right_col - left_col);
            elseif ~isempty(left_vals)
                expected_avg = left_vals;
            elseif ~isempty(right_vals)
                expected_avg = right_vals;
            else
                expected_avg = global_mean;
            end
        end
        corrected(rows_in_col, col) = data(rows_in_col, col) - streak_avg + expected_avg;
    end
end

%% -------------------------------------------------------------------------
function lower_list = thresholdListFromPopulation(L_for_mask, n_thresholds)
% Thresholds at linearly spaced *population* ranks: sort L, take values at rank 1, rank 2, ... (e.g. 10th, 20th, ... point).
% So we get n_thresholds values that span the distribution evenly by count, not by value.
    L_flat = L_for_mask(:);
    L_sorted = sort(L_flat);
    M = numel(L_sorted);
    if M == 0
        lower_list = 0;
        return;
    end
    if n_thresholds >= M
        lower_list = unique(L_sorted, 'stable');
        return;
    end
    % Indices 1, 1 + (M-1)/(N-1), ..., M so we get n_thresholds values (linear in rank)
    indices = round(linspace(1, M, n_thresholds));
    indices = max(1, min(M, indices));
    lower_list = L_sorted(indices);
    lower_list = unique(lower_list, 'stable');
end

%% -------------------------------------------------------------------------
function [var_list, lower_list] = varianceSweep(data, mode, max_streak_width, n_thresholds)
% Sweep thresholds and return variance of corrected image for each (used to pick optimal threshold).
% Uses population-based threshold list (linear in rank, not in value).
    if nargin < 4
        n_thresholds = 100;
    end
    [L_for_mask, ~] = streakCore('laplacian', data, mode, max_streak_width);
    lower_list = thresholdListFromPopulation(L_for_mask, n_thresholds);
    var_list = zeros(size(lower_list));
    for i = 1:numel(lower_list)
        mask = streakCore('mask', L_for_mask, mode, lower_list(i));
        corrected = streakCore('correct', data, mask, 'neighbor_interp');
        var_list(i) = var(corrected(:));
    end
end
