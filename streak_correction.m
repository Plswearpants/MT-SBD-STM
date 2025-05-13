function [corrected_list, var_list, lower_list] = streak_correction(data, mode)
% STREAK_CORRECTION Corrects streaks in image data using Laplacian-based detection
%   [corrected_list, var_list, lower_list] = streak_correction(data, mode)
%
% Inputs:
%   data - 2D matrix containing the image data
%   mode - String specifying correction mode:
%          'both' - Use both left and right neighbors (default)
%          'left' - Use only left neighbors
%
% Outputs:
%   corrected_list - 3D matrix containing corrected images for different thresholds
%   var_list - Vector containing variance values for each threshold
%   lower_list - Vector containing threshold values used
%
% Example:
%   data = Y(:,:,7);
%   [corrected_list, var_list, lower_list] = streak_correction(data);  % Use both neighbors
%   [corrected_list, var_list, lower_list] = streak_correction(data, 'left');  % Use only left neighbors

% Set default mode
if nargin < 2
    mode = 'both';
end

[rows, cols] = size(data);

% Compute Laplacian
L = zeros(size(data));
% Shift left and right
data_left = [zeros(rows,1), data(:,1:end-1)];
data_right = [data(:,2:end), zeros(rows,1)];
% Compute Laplacian using matrix operations
L = data_left + data_right - 2*data;
L_mag = abs(L);

% Generate threshold list
min_mag = min(L_mag(:));
max_mag = max(L_mag(:));
lower_list = linspace(min_mag, max_mag, 100);
var_list = zeros(size(lower_list));

% Initialize corrected_list
corrected_list = zeros(rows, cols, length(lower_list));

% Process each threshold
for i = 1:length(lower_list)
    % Find streaks
    streak_mask = L_mag >= lower_list(i);
    [streak_rows, streak_cols] = find(streak_mask);
    
    % Correct streaks
    corrected = data;
    unique_cols = unique(streak_cols);
    
    for col = unique_cols'
        rows_in_col = streak_rows(streak_cols == col);
        
        % Process all streak points in this column
        if ~isempty(rows_in_col)
            neighbor_cols = [];
            
            % Find left neighbor
            left_col = col - 1;
            if left_col >= 1
                neighbor_cols = [neighbor_cols, left_col];
            end
            
            % Find right neighbor if mode is 'both'
            if strcmpi(mode, 'both')
                right_col = col + 1;
                if right_col <= size(data, 2)
                    neighbor_cols = [neighbor_cols, right_col];
                end
            end
            
            % Get neighbor values
            left_vals = [];
            right_vals = [];
            
            % Process left neighbor if it exists
            if left_col >= 1
                left_vals = mean(data(rows_in_col, left_col));
            end
            
            % Process right neighbor if it exists and mode is 'both'
            if strcmpi(mode, 'both') && col < size(data, 2)
                right_vals = mean(data(rows_in_col, col + 1));
            end
            
            streak_avg = mean(data(rows_in_col, col));
            
            % Calculate expected value based on mode
            if strcmpi(mode, 'both')
                if ~isempty(left_vals) && ~isempty(right_vals)
                    % If both neighbors exist, do linear interpolation
                    left_avg = left_vals;
                    right_avg = right_vals;
                    % Linear interpolation: y = y1 + (x-x1)*(y2-y1)/(x2-x1)
                    expected_avg = left_avg + (col-left_col)*(right_avg-left_avg)/((col+1)-left_col);
                elseif ~isempty(left_vals)
                    % If only left neighbor exists
                    expected_avg = left_vals;
                elseif ~isempty(right_vals)
                    % If only right neighbor exists
                    expected_avg = right_vals;
                else
                    % If no neighbors exist
                    expected_avg = streak_avg;
                end
            else % mode is 'left'
                if ~isempty(left_vals)
                    % If left neighbor exists
                    expected_avg = left_vals;
                else
                    % If no left neighbor exists
                    expected_avg = streak_avg;
                end
            end
            
            corrected(rows_in_col, col) = data(rows_in_col, col) - streak_avg + expected_avg;
        end
    end
    
    % Store results
    corrected_list(:,:,i) = corrected;
    var_list(i) = var(corrected(:));
end

end 