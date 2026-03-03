function [data_out, streak_params] = streakRemovalWorkflow(data_carried, interactive, opts)
%STREAKREMOVALWORKFLOW Single-reference-slice streak detection, global application.
%
%   New minimal workflow:
%   1) User picks one reference slice and tunes threshold interactively.
%   2) The resulting binary streak mask is reused for all slices.
%
%   [data_out, streak_params] = streakRemovalWorkflow(data_carried, interactive)
%   [data_out, streak_params] = streakRemovalWorkflow(data_carried, interactive, opts)
%
%   opts (optional):
%     .streak_params - replay parameters (see output format below)
%
%   streak_params (output) for checkpoint replay:
%     .nSlices          - number of slices
%     .ref_idx          - reference slice used for detection
%     .max_streak_width - detection width on reference slice
%     .mode             - 'valley' | 'plateau' | 'both'
%     .min_value        - final threshold selected on reference slice
%     .streak_mask      - detected reference mask reused for all slices

    if nargin < 3 || ~isstruct(opts)
        opts = struct();
    end

    nSlices = size(data_carried, 3);
    data_out = data_carried;

    if interactive
        ref_idx = round(nSlices / 2);
        max_streak_width = 3;
        mode = 'valley';

        figure; d3gridDisplay(data_carried, 'dynamic');
        fprintf('\nReference slice for streak detection [%d]: ', ref_idx);
        r = input('', 's');
        if ~isempty(strtrim(r))
            v = str2double(r);
            if ~isnan(v) && v >= 1 && v <= nSlices
                ref_idx = round(v);
            end
        end
        close;

        fprintf('max_streak_width [%d]: ', max_streak_width);
        r = input('', 's');
        if ~isempty(strtrim(r))
            v = str2double(r);
            if ~isnan(v) && v >= 1
                max_streak_width = round(v);
            end
        end

        fprintf('mode (valley/plateau/both) [%s]: ', mode);
        r = strtrim(lower(input('', 's')));
        if strcmp(r, 'plateau') || strcmp(r, 'both') || strcmp(r, 'valley')
            mode = r;
        end

        % Interactive detection on one reference slice (UI slider decides final cutoff)
        [~, streak_mask, final_min_val] = removeLocalStreaks( ...
            data_carried, ref_idx, [], max_streak_width, mode, false);

        % Apply same mask to all slices
        h = waitbar(0, 'Applying reference-slice streak mask...');
        for i = 1:nSlices
            [data_out(:,:,i), ~] = removeLocalStreaks(data_carried, i, [], [], mode, true, streak_mask);
            waitbar(i / nSlices, h, sprintf('Processing %d/%d', i, nSlices));
        end
        close(h);

        streak_params = struct();
        streak_params.nSlices = nSlices;
        streak_params.ref_idx = ref_idx;
        streak_params.max_streak_width = max_streak_width;
        streak_params.mode = mode;
        streak_params.min_value = final_min_val;
        streak_params.streak_mask = streak_mask;
        return;
    end

    % Non-interactive replay
    if ~isfield(opts, 'streak_params') || isempty(opts.streak_params)
        error('streakRemovalWorkflow: non-interactive mode requires opts.streak_params.');
    end
    sp = opts.streak_params;
    sp = normalizeReplayParams(sp, nSlices);

    streak_mask = [];
    if isfield(sp, 'streak_mask') && ~isempty(sp.streak_mask)
        streak_mask = logical(sp.streak_mask);
    else
        [~, streak_mask] = removeLocalStreaks( ...
            data_carried, sp.ref_idx, sp.min_value, sp.max_streak_width, sp.mode, true);
    end
    if ~isequal(size(streak_mask), size(data_carried(:,:,1)))
        error('streakRemovalWorkflow: replay streak_mask size mismatch with data slices.');
    end

    for i = 1:nSlices
        [data_out(:,:,i), ~] = removeLocalStreaks(data_carried, i, [], [], sp.mode, true, streak_mask);
    end

    streak_params = sp;
    streak_params.nSlices = nSlices;
    streak_params.streak_mask = streak_mask;
end

function sp = normalizeReplayParams(sp, nSlices)
% Normalize replay params to the new minimal format and keep old checkpoints usable.
    if ~isstruct(sp)
        error('streakRemovalWorkflow: streak_params must be a struct.');
    end

    if isfield(sp, 'ref_idx')
        ref_idx = sp.ref_idx;
    elseif isfield(sp, 'reference_slice')
        ref_idx = sp.reference_slice;
    elseif isfield(sp, 'per_slice') && iscell(sp.per_slice)
        ref_idx = firstRefFromPerSlice(sp.per_slice, nSlices);
    else
        ref_idx = round(nSlices / 2);
    end
    if isempty(ref_idx) || ~isnumeric(ref_idx)
        ref_idx = round(nSlices / 2);
    end
    ref_idx = max(1, min(nSlices, round(ref_idx)));

    max_sw = 3;
    if isfield(sp, 'max_streak_width') && ~isempty(sp.max_streak_width)
        max_sw = round(sp.max_streak_width);
    elseif isfield(sp, 'per_slice') && iscell(sp.per_slice)
        max_sw = getFromPerSlice(sp.per_slice, 'max_streak_width', max_sw);
    end

    mode = 'valley';
    if isfield(sp, 'mode') && ~isempty(sp.mode)
        mode = sp.mode;
    elseif isfield(sp, 'per_slice') && iscell(sp.per_slice)
        mode = getFromPerSlice(sp.per_slice, 'mode', mode);
    end

    min_value = [];
    if isfield(sp, 'min_value') && ~isempty(sp.min_value)
        min_value = sp.min_value;
    elseif isfield(sp, 'remove_factor')
        min_value = [];
    end

    sp.ref_idx = ref_idx;
    sp.max_streak_width = max_sw;
    sp.mode = mode;
    sp.min_value = min_value;
end

function ref_idx = firstRefFromPerSlice(per_slice, nSlices)
    ref_idx = [];
    for i = 1:min(numel(per_slice), nSlices)
        p = per_slice{i};
        if isempty(p) || ~isstruct(p)
            continue;
        end
        if isfield(p, 'reference_slice') && ~isempty(p.reference_slice)
            ref_idx = p.reference_slice;
            return;
        end
        ref_idx = i;
        return;
    end
end

function val = getFromPerSlice(per_slice, field_name, default_val)
    val = default_val;
    for i = 1:numel(per_slice)
        p = per_slice{i};
        if isempty(p) || ~isstruct(p)
            continue;
        end
        if isfield(p, field_name) && ~isempty(p.(field_name))
            val = p.(field_name);
            return;
        end
    end
end
