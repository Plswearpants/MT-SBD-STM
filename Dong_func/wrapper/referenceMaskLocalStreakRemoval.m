function [Y_local_removed, streak_mask, final_min_val] = referenceMaskLocalStreakRemoval( ...
    Y, ref_idx, min_value, max_streak_width, mode, detect_auto, target_slices, show_waitbar)
%REFERENCEMASKLOCALSTREAKREMOVAL Detect on one slice, correct many slices.
%   [Y_local_removed, streak_mask, final_min_val] = referenceMaskLocalStreakRemoval( ...
%       Y, ref_idx, min_value, max_streak_width, mode, detect_auto, target_slices, show_waitbar)
%
%   Workflow:
%   1) Detect streak mask on reference slice ref_idx (interactive or auto).
%   2) Reuse that fixed mask to run local streak correction on target_slices.
%
%   Inputs:
%       Y                - 3D data volume
%       ref_idx          - reference slice index for mask detection
%       min_value        - threshold for removeLocalStreaks on ref slice ([] allowed)
%       max_streak_width - max streak width for detection on ref slice
%       mode             - 'valley' | 'plateau' | 'both'
%       detect_auto      - detection mode for ref slice (false=interactive, true=auto)
%       target_slices    - slices to correct with the fixed mask (default: all slices)
%       show_waitbar     - show progress bar while correcting target_slices (default: false)
%
%   Outputs:
%       Y_local_removed  - corrected 3D volume
%       streak_mask      - mask detected from reference slice
%       final_min_val    - final threshold returned by interactive ref detection

    if nargin < 3
        min_value = [];
    end
    if nargin < 4 || isempty(max_streak_width)
        max_streak_width = 3;
    end
    if nargin < 5 || isempty(mode)
        mode = 'plateau';
    end
    if nargin < 6 || isempty(detect_auto)
        detect_auto = false;
    end
    if nargin < 7 || isempty(target_slices)
        target_slices = 1:size(Y,3);
    end
    if nargin < 8 || isempty(show_waitbar)
        show_waitbar = false;
    end

    nSlices = size(Y, 3);
    if ref_idx < 1 || ref_idx > nSlices
        error('referenceMaskLocalStreakRemoval: ref_idx must be in [1, %d].', nSlices);
    end
    if any(target_slices < 1 | target_slices > nSlices)
        error('referenceMaskLocalStreakRemoval: target_slices must be in [1, %d].', nSlices);
    end
    target_slices = unique(target_slices(:)');

    Y_local_removed = Y;

    % Detect mask on reference slice once
    [Y_local_removed(:,:,ref_idx), streak_mask, final_min_val] = removeLocalStreaks( ...
        Y, ref_idx, min_value, max_streak_width, mode, detect_auto);

    % Apply same mask to all requested slices
    h = [];
    if show_waitbar
        h = waitbar(0, 'Applying reference-slice streak mask...');
    end
    nTarget = numel(target_slices);
    for k = 1:nTarget
        s = target_slices(k);
        if s ~= ref_idx
            [Y_local_removed(:,:,s), ~] = removeLocalStreaks(Y, s, [], [], mode, true, streak_mask);
        end
        if show_waitbar
            waitbar(k / nTarget, h, sprintf('Processing %d/%d (slice %d)', k, nTarget, s));
        end
    end
    if ~isempty(h) && isvalid(h)
        close(h);
    end
end

