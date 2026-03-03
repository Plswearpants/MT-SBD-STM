function [data_out, streak_params] = streakRemovalWorkflow(data_carried, interactive, opts)
%STREAKREMOVALWORKFLOW Per-slice local streak removal: test in a loop, then apply once at End.
%
%   Workflow:
%   1. Maintain a structure: total slice count, each slice -> params (or none).
%   2. Loop: user picks slice list L and params P; manual test on first slice of L;
%      record that slices in L use params P (overwrite per-slice structure for L).
%   3. When user chooses End: apply the final per-slice structure to the volume
%      (each slice with params gets that processing; others unchanged) and output.
%
%   [data_out, streak_params] = streakRemovalWorkflow(data_carried, interactive)
%   [data_out, streak_params] = streakRemovalWorkflow(data_carried, interactive, opts)
%
%   opts (optional):
%     .streak_params  - For replay: struct with .nSlices and .per_slice (see below).
%
%   streak_params (output) for checkpoint replay:
%     .nSlices   - Number of slices in the volume.
%     .per_slice - cell(1, nSlices). per_slice{s} = [] means no processing;
%                  per_slice{s} = struct(.max_streak_width, .mode, .remove_factor,
%                  .interpolate_factor) means apply that to slice s.
%   Replay: caller passes opts.streak_params; this function applies per_slice to data (no UI).

    if nargin < 3
        opts = struct();
    end
    if ~isstruct(opts)
        opts = struct();
    end

    nSlices = size(data_carried, 3);

    % ----- Replay: apply per_slice structure to data (no UI) -----
    if ~interactive
        if ~isfield(opts, 'streak_params') || isempty(opts.streak_params)
            error('streakRemovalWorkflow: non-interactive mode requires opts.streak_params.');
        end
        sp = opts.streak_params;
        nSlices_data = size(data_carried, 3);
        if isfield(sp, 'per_slice') && iscell(sp.per_slice) && numel(sp.per_slice) >= nSlices_data
            % New format: apply per_slice
            data_out = applyPerSliceStructure(data_carried, sp.per_slice, nSlices_data);
        else
            % Backward compat: single .slices + single params -> build per_slice
            sp = legacyStreakParamsToPerSlice(sp, nSlices_data);
            data_out = applyPerSliceStructure(data_carried, sp.per_slice, nSlices_data);
        end
        streak_params = sp;
        if ~isfield(streak_params, 'nSlices')
            streak_params.nSlices = nSlices_data;
        end
        return;
    end

    % ----- Interactive: per-slice structure, overwritten each loop -----
    data_before_streak = data_carried;
    per_slice = cell(1, nSlices);
    for i = 1:nSlices
        per_slice{i} = [];
    end

    % Default parameters (exposed each loop)
    max_streak_width = 3;
    mode = 'plateau';
    remove_factor = 1;
    interpolate_factor = 0.8;

    done = false;
    while ~done
        figure; d3gridDisplay(data_before_streak, 'dynamic');
        % --- Slice list for this run ---
        fprintf('\nEnter slice list for this run (e.g. 1:50 or [1 5 10 20]): ');
        inp = input('', 's');
        if isempty(strtrim(inp))
            slices = 1:min(50, nSlices);
        else
            try
                slices = eval(inp);
                slices = slices(:)';
            catch
                error('streakRemovalWorkflow: invalid slice list. Use e.g. 1:50 or [1 5 10].');
            end
        end
        if isempty(slices) || any(slices < 1 | slices > nSlices)
            error('streakRemovalWorkflow: slices must be in range [1, %d].', nSlices);
        end
        close();
        % --- Parameters for this run (overwrites per_slice for these slices) ---
        fprintf('  max_streak_width [%d]: ', max_streak_width);
        r = input('', 's');
        if ~isempty(strtrim(r))
            v = str2double(r);
            if ~isnan(v) && v >= 1, max_streak_width = round(v); end
        end
        fprintf('  mode (plateau/valley) [%s]: ', mode);
        r = input('', 's');
        if ~isempty(strtrim(r))
            if strcmpi(strtrim(r), 'valley')
                mode = 'valley';
            else
                mode = 'plateau';
            end
        end
        fprintf('  remove_factor [%.2f]: ', remove_factor);
        r = input('', 's');
        if ~isempty(strtrim(r))
            v = str2double(r);
            if ~isnan(v), remove_factor = v; end
        end
        fprintf('  interpolate_factor [%.2f]: ', interpolate_factor);
        r = input('', 's');
        if ~isempty(strtrim(r))
            v = str2double(r);
            if ~isnan(v), interpolate_factor = v; end
        end

        P = struct();
        P.max_streak_width = max_streak_width;
        P.mode = mode;
        P.remove_factor = remove_factor;
        P.interpolate_factor = interpolate_factor;

        % --- Manual test on first slice of L: sliders set initial position from factors; on Done we read back effective factors ---
        s_first = slices(1);
        [~, var_list, low_list] = streak_correction(data_before_streak(:,:,s_first), max_streak_width, mode);
        figure('Name', 'Variance vs threshold (first slice)');
        plot(low_list, var_list);
        xlabel('Threshold'); ylabel('Variance');
        title(sprintf('Slice %d - adjust sliders then click Done', s_first));
        [~, min_idx] = min(var_list);
        min_low = low_list(min_idx);
        if min_low == 0
            min_low = eps;
        end
        % Remove step: interactive; initial threshold = remove_factor*min_low; returns final threshold on Done
        [Y_after_remove, ~, final_min_val_remove] = removeLocalStreaks(data_before_streak, s_first, remove_factor*min_low, ...
            max_streak_width, mode, false);
        % Interpolate step: interactive; initial threshold = interpolate_factor*min_low; returns final threshold on Done
        [Y_first, ~, ~, final_min_val_interp] = interpolateLocalStreaks(Y_after_remove, 1, interpolate_factor*min_low, [], false);
        % Compute effective factors from slider positions (feedback for next run and for recording)
        if ~isempty(final_min_val_remove)
            remove_factor = final_min_val_remove / min_low;
            fprintf('  [from slider] remove_factor = %.4f\n', remove_factor);
        end
        if ~isempty(final_min_val_interp)
            interpolate_factor = final_min_val_interp / min_low;
            fprintf('  [from slider] interpolate_factor = %.4f\n', interpolate_factor);
        end
        % Use these effective factors for P (what user actually chose)
        P.remove_factor = remove_factor;
        P.interpolate_factor = interpolate_factor;
        % Show result
        figure('Name', sprintf('Result slice %d (params for this run)', s_first));
        imagesc(Y_first); axis image; colormap(gca, parula); colorbar;

        % --- Overwrite per_slice for slices in L with P (P now has slider-derived factors) ---
        for s = slices
            per_slice{s} = P;
        end
        fprintf('Recorded params for slices %s (remove_factor=%.4f, interpolate_factor=%.4f). Total slices with params: %d.\n', ...
            mat2str(slices), P.remove_factor, P.interpolate_factor, sum(cellfun(@(x) ~isempty(x), per_slice)));

        % --- Repeat or end ---
        fprintf('Repeat (test more slices / different params) or End? [r/e]: ');
        reply = input('', 's');
        if ~strcmpi(strtrim(reply), 'r') && ~strcmpi(strtrim(reply), 'repeat')
            done = true;
        end
    end

    % --- Apply final per_slice structure to full volume and output ---
    data_out = applyPerSliceStructure(data_before_streak, per_slice, nSlices);

    streak_params = struct();
    streak_params.nSlices = nSlices;
    streak_params.per_slice = per_slice;
end

function data_out = applyPerSliceStructure(data_in, per_slice, nSlices)
% Apply per-slice params: for each s with non-empty per_slice{s}, run remove+interpolate.
    data_out = data_in;
    for s = 1:nSlices
        if s > numel(per_slice) || isempty(per_slice{s})
            continue;
        end
        P = per_slice{s};
        if ~isstruct(P) || ~isfield(P, 'max_streak_width')
            continue;
        end
        max_sw = P.max_streak_width;
        mode = 'plateau';
        if isfield(P, 'mode'), mode = P.mode; end
        remove_fac = 1;
        if isfield(P, 'remove_factor'), remove_fac = P.remove_factor; end
        interp_fac = 1;
        if isfield(P, 'interpolate_factor'), interp_fac = P.interpolate_factor; end
        [~, var_list, low_list] = streak_correction(data_out(:,:,s), max_sw, mode);
        [~, min_idx] = min(var_list);
        min_low = low_list(min_idx);
        [data_out(:,:,s), ~] = removeLocalStreaks(data_out, s, remove_fac*min_low, max_sw, mode, true);
        [data_out(:,:,s), ~] = interpolateLocalStreaks(data_out(:,:,s), 1, interp_fac*min_low, [], true);
        fprintf('Applied params to slice %d\n', s);
    end
end

function sp = legacyStreakParamsToPerSlice(sp, nSlices)
% Convert old format (.slices + single params) to .per_slice.
    per_slice = cell(1, nSlices);
    for i = 1:nSlices
        per_slice{i} = [];
    end
    if ~isfield(sp, 'slices') || isempty(sp.slices)
        sp.per_slice = per_slice;
        sp.nSlices = nSlices;
        return;
    end
    P = struct();
    P.max_streak_width = 3;
    if isfield(sp, 'max_streak_width'), P.max_streak_width = sp.max_streak_width; end
    P.mode = 'plateau';
    if isfield(sp, 'mode'), P.mode = sp.mode; end
    P.remove_factor = 1;
    if isfield(sp, 'remove_factor'), P.remove_factor = sp.remove_factor; end
    P.interpolate_factor = 1;
    if isfield(sp, 'interpolate_factor'), P.interpolate_factor = sp.interpolate_factor; end
    for s = sp.slices
        if s >= 1 && s <= nSlices
            per_slice{s} = P;
        end
    end
    sp.per_slice = per_slice;
    sp.nSlices = nSlices;
end
