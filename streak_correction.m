function [corrected_list, var_list, lower_list] = streak_correction(data, max_streak_width, mode)
% STREAK_CORRECTION Corrects streaks in image data using multi-scale Laplacian-based detection
%   [corrected_list, var_list, lower_list] = streak_correction(data, max_streak_width, mode)
%
% Inputs:
%   data - 2D matrix containing the image data
%   max_streak_width - Maximum streak width to detect (default: 3)
%   mode - 'valley', 'plateau', or 'both' (default: 'both')
%
% Outputs:
%   corrected_list - 3D matrix containing corrected images for different thresholds
%   var_list - Vector containing variance values for each threshold
%   lower_list - Vector containing threshold values used
%
% Example:
%   data = Y(:,:,7);
%   [corrected_list, var_list, lower_list] = streak_correction(data);  % Default max width 3, both
%   [corrected_list, var_list, lower_list] = streak_correction(data, 5, 'valley');  % Max width 5, valley mode

% Set default parameters
if nargin < 2
    max_streak_width = 3;
end
if nargin < 3
    mode = 'both';
end

[rows, cols] = size(data);

% Compute global mean for fallback
global_mean = mean(data(:));

% Compute multi-scale Laplacians
L_combined = zeros(size(data));
L_mag_combined = zeros(size(data));

for n = 1:max_streak_width
    % Shift left and right by n pixels
    data_left_n = [zeros(rows, n), data(:, 1:end-n)];
    data_right_n = [data(:, n+1:end), zeros(rows, n)];
    
    % Compute Laplacian at distance n
    L_n = data_left_n + data_right_n - 2*data;
    
    % Combine with previous scales (take maximum magnitude)
    if strcmpi(mode, 'valley') || strcmpi(mode, 'plateau')
        L_combined = max(L_combined, L_n); % for threshold range
    else
        L_mag_n = abs(L_n);
        L_mag_combined = max(L_mag_combined, L_mag_n);
    end
end

% Generate threshold list
if strcmpi(mode, 'valley')
    L_for_mask = L_combined;
    min_thr = 0;
    max_thr = max(L_for_mask(:));
    lower_list = linspace(min_thr, max_thr, 300);
elseif strcmpi(mode, 'plateau')
    L_for_mask = L_combined;
    min_thr = min(L_for_mask(:));
    max_thr = 0;
    lower_list = linspace(min_thr, max_thr, 300);
else % both
    L_for_mask = L_mag_combined;
    min_thr = min(L_for_mask(:));
    max_thr = prctile(L_for_mask(:), 99.5);
    lower_list = linspace(min_thr, max_thr, 300);
end
var_list = zeros(size(lower_list));

% Initialize corrected_list
corrected_list = zeros(rows, cols, length(lower_list));

% Process each threshold
for i = 1:length(lower_list)
    if strcmpi(mode, 'valley')
        streak_mask = (L_for_mask > lower_list(i));
    elseif strcmpi(mode, 'plateau')
        streak_mask = (L_for_mask < lower_list(i));
    else
        streak_mask = (L_for_mask >= lower_list(i));
    end
    [streak_rows, streak_cols] = find(streak_mask);
    
    % Correct streaks
    corrected = data;
    unique_cols = unique(streak_cols);
    
    for col = unique_cols'
        rows_in_col = streak_rows(streak_cols == col);
        
        % Process all streak points in this column
        if ~isempty(rows_in_col)
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
            streak_avg = mean(data(rows_in_col, col));
            % Calculate expected value
            if ~isempty(left_vals) && ~isempty(right_vals)
                % If both neighbors exist, do linear interpolation
                left_avg = left_vals;
                right_avg = right_vals;
                expected_avg = left_avg + (col-left_col)*(right_avg-left_avg)/(right_col-left_col);
            elseif ~isempty(left_vals)
                expected_avg = left_vals;
            elseif ~isempty(right_vals)
                expected_avg = right_vals;
            else
                expected_avg = global_mean;
            end
            corrected(rows_in_col, col) = data(rows_in_col, col) - streak_avg + expected_avg;
        end
    end
    
    % Store results
    corrected_list(:,:,i) = corrected;
    var_list(i) = var(corrected(:));
end

end 