function [corrected_list, var_list, lower_list] = streak_correction(data, max_streak_width, mode)
% STREAK_CORRECTION Corrects streaks in image data using multi-scale Laplacian-based detection
%   [corrected_list, var_list, lower_list] = streak_correction(data, max_streak_width, mode)
%
%   Threshold list is built from the *population* (sorted Laplacian values): we take
%   thresholds at linearly spaced ranks (e.g. 10th, 20th, ... point), not linear value
%   spacing, so resolution follows the data distribution.
%
% Inputs:
%   data - 2D matrix
%   max_streak_width - default 3
%   mode - 'valley', 'plateau', or 'both'
%
% Outputs:
%   corrected_list - 3D array of corrected images per threshold
%   var_list - variance of each corrected image
%   lower_list - threshold values (pick min_low = lower_list(argmin(var_list)) for best threshold)

if nargin < 2
    max_streak_width = 3;
end
if nargin < 3
    mode = 'both';
end

[rows, cols] = size(data);
[L_for_mask, ~] = streakCore('laplacian', data, mode, max_streak_width);
lower_list = streakCore('threshold_list', L_for_mask, 100);
var_list = zeros(size(lower_list));
corrected_list = zeros(rows, cols, numel(lower_list));

for i = 1:numel(lower_list)
    mask = streakCore('mask', L_for_mask, mode, lower_list(i));
    corrected = streakCore('correct', data, mask, 'neighbor_interp');
    corrected_list(:,:,i) = corrected;
    var_list(i) = var(corrected(:));
end
end
