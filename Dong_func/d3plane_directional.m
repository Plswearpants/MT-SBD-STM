function [data_plane, mask] = d3plane_directional(data, direction, varargin)
% this function correct the data to make it plane on a specific direction, 
% via manual drawing linesegments along that direction. The logic is that points on the 
% line segments should have the same value. 
%
% Input:
%   data: 3D data array
%   direction: 'horizontal' or 'vertical'
%   'LineWidth': (optional) Width of the line segments in pixels (default: 1)
%
% Output:
%   data_plane: 3D data array that is planed on the given direction 

% Parse optional inputs
p = inputParser;
addParameter(p, 'LineWidth', 1, @(x) isnumeric(x) && x > 0);
parse(p, varargin{:});
line_width = p.Results.LineWidth;

% Get mask from line segments
mask = linesegmentEdge2Edge(data, direction, 'LineWidth', line_width);

% Create NaN mask for non-selected regions
data_offset = data;
data_offset(~logical(mask)) = NaN;

switch direction
    case 'horizontal'
        % For each row, calculate the mean value along the masked region, here the size of row_mean is [x,1,slice]
        row_means = mean(data_offset, 1, 'omitnan');
        
        % Debug: Check for NaNs in row_means
        if any(isnan(row_means(:)))
            fprintf('Warning: NaNs found in row_means\n');
            fprintf('Number of NaNs: %d\n', sum(isnan(row_means(:))));
            fprintf('Size of row_means: %s\n', mat2str(size(row_means)));
            fprintf('Min value: %f\n', min(row_means(:), [], 'omitnan'));
            fprintf('Max value: %f\n', max(row_means(:), [], 'omitnan'));
        end
        
        % apply 1D gaussian smoothing along dimension 1 to the row_means
        row_means = smoothdata(row_means, 1, 'gaussian', 5, 'omitnan');
        
        % Calculate the reference value (e.g., mean of all row means)
        ref_value = mean(row_means, 'omitnan');
        
        % Calculate offsets for each row
        row_offsets = row_means - ref_value;

        % repeat the row_means to the size of data
        row_means = repmat(row_means, [size(data, 1), 1, 1]);
        
        % Apply offsets to all slices
        data_plane = data - row_means;
        
    case 'vertical'
        % For each column, calculate the mean value along the masked region
        col_means = mean(data_offset, 2, 'omitnan');
        
        % Debug: Check for NaNs in col_means
        if any(isnan(col_means(:)))
            fprintf('Warning: NaNs found in col_means\n');
            fprintf('Number of NaNs: %d\n', sum(isnan(col_means(:))));
            fprintf('Size of col_means: %s\n', mat2str(size(col_means)));
            fprintf('Min value: %f\n', min(col_means(:), [], 'omitnan'));
            fprintf('Max value: %f\n', max(col_means(:), [], 'omitnan'));
        end
        
        % apply 1D gaussian smoothing along dimension 2 to the col_means
        col_means = smoothdata(col_means, 2, 'gaussian', 5, 'omitnan');
        
        % Calculate the reference value (e.g., mean of all column means)
        ref_value = mean(col_means, 'omitnan');
        
        % Calculate offsets for each column
        col_offsets = col_means - ref_value;
        
        % repeat the col_means to the size of data
        col_means = repmat(col_means, [1, size(data, 2), 1]);
        
        % Apply offsets to all slices
        data_plane = data - col_means;
end

% Display results
figure;
subplot(1,3,1);
imagesc(data(:,:,1));
title('Original Data');
colormap(gray);
axis square;

subplot(1,3,2);
imagesc(mask);
title(sprintf('Mask (Line Width: %d pixels)', line_width));
colormap(gray);
axis square;

subplot(1,3,3);
imagesc(data_plane(:,:,1));
title('Corrected Data');
colormap(gray);
axis square;

end
