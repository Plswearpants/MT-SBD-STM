function [mask_ref, chosen_threshold, L_for_mask] = singleSliceStreakMask(Y, ref_idx, mode, max_streak_width, n_thresholds)
%SINGLESLICESTREAKMASK Compute streak mask from a single reference slice.
%   [mask_ref, chosen_threshold, L_for_mask] = singleSliceStreakMask(Y, ref_idx, mode, max_streak_width, n_thresholds)
%
%   Inputs:
%       Y               - 3D volume (rows x cols x slices)
%       ref_idx         - reference slice index to run detection on
%       mode            - 'valley', 'plateau', or 'both' (default: 'plateau')
%       max_streak_width - maximum streak width (default: 3)
%       n_thresholds    - number of thresholds to sweep for variance-based choice (default: 100)
%
%   Outputs:
%       mask_ref         - binary streak mask from the reference slice (logical matrix)
%       chosen_threshold - scalar threshold used to build mask_ref
%       L_for_mask       - Laplacian image used for mask construction (for diagnostics)
%
%   The threshold is selected by sweeping thresholds with streakCore('sweep', ...)
%   and picking the one that minimizes the variance of the corrected image.

    if nargin < 3 || isempty(mode)
        mode = 'plateau';
    end
    if nargin < 4 || isempty(max_streak_width)
        max_streak_width = 3;
    end
    if nargin < 5 || isempty(n_thresholds)
        n_thresholds = 100;
    end

    if ndims(Y) ~= 3
        error('singleSliceStreakMask: Y must be a 3D volume (rows x cols x slices).');
    end
    nSlices = size(Y, 3);
    if ref_idx < 1 || ref_idx > nSlices
        error('singleSliceStreakMask: ref_idx must be in [1, %d].', nSlices);
    end

    data_ref = Y(:,:,ref_idx);

    % Compute Laplacian field for the reference slice
    [L_for_mask, ~] = streakCore('laplacian', data_ref, mode, max_streak_width);

    % Sweep thresholds using shared core logic and pick the one that minimizes variance
    [var_list, lower_list] = streakCore('sweep', data_ref, mode, max_streak_width, n_thresholds);
    [~, min_idx] = min(var_list);
    chosen_threshold = lower_list(min_idx);

    % Final binary mask from the chosen threshold (includes vertical continuity filter)
    mask_ref = streakCore('mask', L_for_mask, mode, chosen_threshold);
end
