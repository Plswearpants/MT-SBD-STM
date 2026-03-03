function Y_corr = applyStreakMaskAllSlices(Y, mask_ref, mode, method)
%APPLYSTREAKMASKALLSLICES Apply a fixed streak mask to all slices of a volume.
%   Y_corr = applyStreakMaskAllSlices(Y, mask_ref, mode, method)
%
%   Inputs:
%       Y        - 3D volume (rows x cols x slices)
%       mask_ref - binary streak mask (logical or numeric) from a reference slice
%       mode     - 'valley', 'plateau', or 'both' (used to choose default correction method)
%       method   - correction method string for streakCore('correct', ...):
%                  'neighbor_interp' | 'background_mean' | 'horizontal_avg'
%                  (optional; if omitted, picks based on mode, matching existing logic)
%
%   Output:
%       Y_corr   - corrected 3D volume, same size as Y

    if ndims(Y) ~= 3
        error('applyStreakMaskAllSlices: Y must be a 3D volume (rows x cols x slices).');
    end

    [rows, cols, nSlices] = size(Y);
    if ~isequal(size(mask_ref), [rows, cols])
        error('applyStreakMaskAllSlices: mask_ref size must match first two dimensions of Y.');
    end

    if nargin < 3 || isempty(mode)
        mode = 'plateau';
    end

    if nargin < 4 || isempty(method)
        % Match the existing behavior in streakCore('correct', ...)
        if strcmpi(mode, 'valley') || strcmpi(mode, 'plateau')
            method = 'background_mean';
        else
            method = 'neighbor_interp';
        end
    end

    mask_ref = logical(mask_ref);
    Y_corr = Y;

    for s = 1:nSlices
        data_s = Y(:,:,s);
        corrected_s = streakCore('correct', data_s, mask_ref, method);
        Y_corr(:,:,s) = corrected_s;
    end
end

function Y_corr = applyStreakMaskAllSlices(Y, mask_ref, mode, method)
%APPLYSTREAKMASKALLSLICES Apply a fixed streak mask to all slices of a volume.
%   Y_corr = applyStreakMaskAllSlices(Y, mask_ref, mode, method)
%
%   Inputs:
%       Y        - 3D volume (rows x cols x slices)
%       mask_ref - binary streak mask (logical or numeric) from a reference slice
%       mode     - 'valley', 'plateau', or 'both' (used to choose default correction method)
%       method   - correction method string for streakCore('correct', ...):
%                  'neighbor_interp' | 'background_mean' | 'horizontal_avg'
%                  (optional; if omitted, picks based on mode, matching existing logic)
%
%   Output:
%       Y_corr   - corrected 3D volume, same size as Y

    if ndims(Y) ~= 3
        error('applyStreakMaskAllSlices: Y must be a 3D volume (rows x cols x slices).');
    end

    [rows, cols, nSlices] = size(Y);
    if ~isequal(size(mask_ref), [rows, cols])
        error('applyStreakMaskAllSlices: mask_ref size must match first two dimensions of Y.');
    end

    if nargin < 3 || isempty(mode)
        mode = 'plateau';
    end

    if nargin < 4 || isempty(method)
        % Match the existing behavior in streakCore('correct', ...)
        if strcmpi(mode, 'valley') || strcmpi(mode, 'plateau')
            method = 'background_mean';
        else
            method = 'neighbor_interp';
        end
    end

    mask_ref = logical(mask_ref);
    Y_corr = Y;

    for s = 1:nSlices
        data_s = Y(:,:,s);
        corrected_s = streakCore('correct', data_s, mask_ref, method);
        Y_corr(:,:,s) = corrected_s;
    end
end

