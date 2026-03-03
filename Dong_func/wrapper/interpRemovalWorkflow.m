function [data_out, interp_params] = interpRemovalWorkflow(data_in, interactive, opts)
%INTERPREMOVALWORKFLOW Single-reference-slice interpolation of streak pixels.
%
%   Workflow:
%   1) User picks one reference slice and tunes interpolation threshold interactively
%      via interpolateLocalStreaks (UI).
%   2) The resulting interpolation mask (pixels to be in-painted) is reused for all
%      slices: those pixels are set to NaN and filled by natural-neighbor interpolation.
%
%   [data_out, interp_params] = interpRemovalWorkflow(data_in, interactive)
%   [data_out, interp_params] = interpRemovalWorkflow(data_in, interactive, opts)
%
%   Inputs:
%       data_in     - 3D volume (rows x cols x slices)
%       interactive - true for UI-based reference selection, false for replay
%       opts        - optional struct for replay:
%                     .interp_params (see below)
%
%   Outputs:
%       data_out     - interpolated 3D volume
%       interp_params.nSlices      - number of slices
%       interp_params.ref_idx      - reference slice index
%       interp_params.min_value    - final threshold selected on reference slice
%       interp_params.interp_mask  - logical mask of pixels interpolated on ref slice
%       interp_params.streak_indices - [N x 2] row/col indices on ref slice

    if nargin < 3 || ~isstruct(opts)
        opts = struct();
    end

    if ndims(data_in) ~= 3
        error('interpRemovalWorkflow: data_in must be 3D (rows x cols x slices).');
    end
    [rows, cols, nSlices] = size(data_in);

    data_out = data_in;

    if interactive
        ref_idx = round(nSlices / 2);

        figure; d3gridDisplay(data_in, 'dynamic');
        fprintf('\nReference slice for interpolation [%d]: ', ref_idx);
        r = input('', 's');
        if ~isempty(strtrim(r))
            v = str2double(r);
            if ~isnan(v) && v >= 1 && v <= nSlices
                ref_idx = round(v);
            end
        end
        close;

        min_value = [];
        provided_streak_indices = [];
        auto = false;
        [corrected_ref, interp_mask, streak_indices, final_min_val] = interpolateLocalStreaks( ...
            data_in, ref_idx, min_value, provided_streak_indices, auto);

        % Apply same interpolation pattern to all slices using natural-neighbor fill
        mask = logical(interp_mask);
        [X, Y] = meshgrid(1:cols, 1:rows);

        h = waitbar(0, 'Interpolating streak pixels across slices...');
        for k = 1:nSlices
            slice = double(data_in(:,:,k));

            % Mark interpolation pixels as NaN
            slice(mask) = NaN;

            % Natural-neighbor interpolation on non-NaN pixels
            known = ~isnan(slice);
            F = scatteredInterpolant(X(known), Y(known), slice(known), 'natural', 'nearest');
            slice(mask) = F(X(mask), Y(mask));

            data_out(:,:,k) = slice;
            waitbar(k/nSlices, h, sprintf('Processing %d/%d', k, nSlices));
        end
        close(h);

        interp_params = struct();
        interp_params.nSlices = nSlices;
        interp_params.ref_idx = ref_idx;
        interp_params.min_value = final_min_val;
        interp_params.interp_mask = mask;
        interp_params.streak_indices = streak_indices;
        return;
    end

    % Non-interactive replay
    if ~isfield(opts, 'interp_params') || isempty(opts.interp_params)
        error('interpRemovalWorkflow: non-interactive mode requires opts.interp_params.');
    end
    ip = opts.interp_params;
    ip = normalizeInterpReplayParams(ip, nSlices, [rows, cols]);

    mask = logical(ip.interp_mask);
    [X, Y] = meshgrid(1:cols, 1:rows);

    for k = 1:nSlices
        slice = double(data_in(:,:,k));
        slice(mask) = NaN;
        known = ~isnan(slice);
        F = scatteredInterpolant(X(known), Y(known), slice(known), 'natural', 'nearest');
        slice(mask) = F(X(mask), Y(mask));
        data_out(:,:,k) = slice;
    end

    interp_params = ip;
    interp_params.nSlices = nSlices;
end

function ip = normalizeInterpReplayParams(ip, nSlices, spatialSize)
% Normalize interpolation replay parameters to the new minimal format.

    if ~isstruct(ip)
        error('interpRemovalWorkflow: interp_params must be a struct.');
    end

    if isfield(ip, 'ref_idx') && ~isempty(ip.ref_idx)
        ref_idx = ip.ref_idx;
    else
        ref_idx = round(nSlices / 2);
    end
    ref_idx = max(1, min(nSlices, round(ref_idx)));

    min_value = [];
    if isfield(ip, 'min_value') && ~isempty(ip.min_value)
        min_value = ip.min_value;
    end

    if isfield(ip, 'interp_mask') && ~isempty(ip.interp_mask)
        mask = logical(ip.interp_mask);
        if ~isequal(size(mask), spatialSize)
            error('interpRemovalWorkflow: replay interp_mask size mismatch with data slices.');
        end
    else
        error('interpRemovalWorkflow: replay currently requires interp_mask to be stored.');
    end

    ip.ref_idx = ref_idx;
    ip.min_value = min_value;
    ip.interp_mask = mask;
end

