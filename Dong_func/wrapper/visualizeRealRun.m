function [log, data, params, meta, cfg] = visualizeRealRun(log, data, params, meta, cfg)
%VISUALIZEREALRUN Basic visualization of real-data block-run results.
%
%   [log, data, params, meta, cfg] = visualizeRealRun(log, data, params, meta, cfg)
%
%   This is a lightweight wrapper intended for VR01R in run_real_data.m.
%   It reuses existing helpers (visualizeRealResult, d3gridDisplay,
%   qpiCalculate) to provide a quick look at selected slices.

    if ~isfield(data, "real") || ~isfield(data.real, "blockRun")
        warning('visualizeRealRun: block-run results not found; nothing to visualize.');
        return;
    end

    Y = data.real.Y;
    Aout_ALL   = data.real.blockRun.Aout_ALL;
    Xout_ALL   = data.real.blockRun.Xout_ALL;
    bout_ALL   = data.real.blockRun.bout_ALL;
    ALL_extras = data.real.blockRun.ALL_extras;

    [num_slices, num_kernels] = size(bout_ALL);

    % Convert Aout_ALL to per-slice cell format (as in legacy script)
    Aout_ALL_cell = cell(num_slices, num_kernels);
    for s = 1:num_slices
        for k = 1:num_kernels
            Aout_ALL_cell{s,k} = Aout_ALL{k}(:,:,s);
        end
    end

    % Choose slices to visualize
    if isfield(params, "visualize") && isfield(params.visualize, "all_slices") ...
            && ~params.visualize.all_slices
        slices_to_viz = params.refSlice.ref_slice;
    else
        slices_to_viz = 1:num_slices;
    end

    for s = slices_to_viz
        fprintf('Visualizing slice %d/%d...\n', s, num_slices);
        pp = struct();
        pp.phase1 = struct();
        if isfield(ALL_extras, "phase1")
            pp.phase1.residuals        = ALL_extras.phase1.residuals(:,:,s,:);
            pp.phase1.quality_metrics  = ALL_extras.phase1.quality_metrics;
        end

        visualizeRealResult(Y(:,:,s), Aout_ALL_cell(s,:), Xout_ALL, bout_ALL(s,:), pp);
    end

end

