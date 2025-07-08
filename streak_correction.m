function [corrected_list, var_list, lower_list] = streak_correction(data, max_streak_width)
% STREAK_CORRECTION Corrects streaks in image data using multi-scale Laplacian-based detection
%   [corrected_list, var_list, lower_list] = streak_correction(data, max_streak_width)
%
% Inputs:
%   data - 2D matrix containing the image data
%   max_streak_width - Maximum streak width to detect (default: 3)
%
% Outputs:
%   corrected_list - 3D matrix containing corrected images for different thresholds
%   var_list - Vector containing variance values for each threshold
%   lower_list - Vector containing threshold values used
%
% Example:
%   data = Y(:,:,7);
%   [corrected_list, var_list, lower_list] = streak_correction(data);  % Default max width 3
%   [corrected_list, var_list, lower_list] = streak_correction(data, 5);  % Max width 5

% Set default parameters
if nargin < 2
    max_streak_width = 3;
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
    L_mag_n = abs(L_n);
    
    % Combine with previous scales (take maximum magnitude)
    L_mag_combined = max(L_mag_combined, L_mag_n);
end

% Generate threshold list
min_mag = min(L_mag_combined(:));
max_mag = prctile(L_mag_combined(:), 99.5);
lower_list = linspace(min_mag, max_mag, 300);
var_list = zeros(size(lower_list));

% Initialize corrected_list
corrected_list = zeros(rows, cols, length(lower_list));

% Process each threshold
for i = 1:length(lower_list)
    % Find streaks using combined magnitude
    streak_mask = L_mag_combined >= lower_list(i);
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