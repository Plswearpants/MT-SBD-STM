function [log, data, params, meta, cfg] = visualizeRealRun(log, data, params, meta, cfg)
%VISUALIZEREALRUN Visualization of real-data block-run results.
%
%   [log, data, params, meta, cfg] = visualizeRealRun(log, data, params, meta, cfg)
%
%   Intended for VR01R in run_real_data.m.
%   It produces:
%       (1) Kernel movies: for each kernel, d3gridDisplay across slices.
%       (2) Reconstruction movie: conv(Aout, Xout) + bias for all slices.
%
    if ~isfield(data, "real") || ~isfield(data.real, "blockRun")
        warning('visualizeRealRun: block-run results not found; nothing to visualize.');
        return;
    end

    Y = data.real.Y;
    Aout_ALL   = data.real.blockRun.Aout_ALL;
    Xout_ALL   = data.real.blockRun.Xout_ALL;
    bout_ALL   = data.real.blockRun.bout_ALL;

    [num_slices_run, num_kernels] = size(bout_ALL);

    % If block run used a subset, map results back to original slice indices.
    if isfield(data.real.blockRun, "slice_indices") && ~isempty(data.real.blockRun.slice_indices)
        slice_indices = data.real.blockRun.slice_indices(:).';
    else
        slice_indices = 1:num_slices_run;
    end

    % Convert Aout_ALL to per-slice cell format
    Aout_ALL_cell = cell(num_slices_run, num_kernels);
    for s = 1:num_slices_run
        for k = 1:num_kernels
            Aout_ALL_cell{s,k} = Aout_ALL{k}(:,:,s);
        end
    end

    %----------------------------------------------------------------------
    % 1) Kernel movies
    %----------------------------------------------------------------------
    for k = 1:num_kernels
        fprintf('Kernel movie for kernel %d/%d...\n', k, num_kernels);
        figure;
        d3gridDisplay(Aout_ALL{k}, 'dynamic');
        title(sprintf('Kernel %d across slices', k));
    end

    %----------------------------------------------------------------------
    % 2) Reconstruction movie for all slices
    %    Y_rec(:,:,s_full) = sum_k convfft2(Aout(s_full,k), Xout(:,:,k)) + bias
    %----------------------------------------------------------------------
    Y_rec_full = zeros([size(Y,[1,2]),num_slices_run]);
    for idx = 1:num_slices_run
        for k = 1:num_kernels
            Y_rec_full(:,:,idx) = Y_rec_full(:,:,idx) + convfft2(Aout_ALL_cell{idx,k}, Xout_ALL(:,:,k)) + bout_ALL(idx,k);
        end
    end

    figure;
    d3gridDisplay(Y_rec_full, 'dynamic');
    title('Reconstructed observation (all slices)');

end

