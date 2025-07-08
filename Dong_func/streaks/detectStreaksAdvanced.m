function [streak_mask, streak_components, visualization_data] = detectStreaksAdvanced(data, varargin)
%DETECTSTREAKSADVANCED Advanced streak detection using edge-pair analysis
%   [streak_mask, streak_components, visualization_data] = detectStreaksAdvanced(data, varargin)
%
% Inputs:
%   data - 2D matrix containing the image data
%   'Threshold' - Laplacian threshold (default: auto-determined)
%   'MaxStreakWidth' - Maximum width between edge pairs (default: 3)
%   'ShowVisualization' - Show interactive visualization (default: true)
%
% Outputs:
%   streak_mask - Binary mask indicating detected streaks
%   streak_components - Structure array with component information
%   visualization_data - Structure containing visualization data
%
% Example:
%   [mask, components, viz] = detectStreaksAdvanced(data);
%   [mask, components, viz] = detectStreaksAdvanced(data, 'Threshold', 0.5, 'MaxStreakWidth', 5);

% Parse optional inputs
p = inputParser;
addParameter(p, 'Threshold', [], @isnumeric);
addParameter(p, 'MaxStreakWidth', 3, @isnumeric);
addParameter(p, 'ShowVisualization', true, @islogical);
parse(p, varargin{:});

threshold = p.Results.Threshold;
max_width = p.Results.MaxStreakWidth;
show_viz = p.Results.ShowVisualization;

[rows, cols] = size(data);

%% Step 1: Compute Laplacian
fprintf('Computing Laplacian...\n');

% Compute Laplacian
L = zeros(size(data));
% Shift left and right
data_left = [zeros(rows,1), data(:,1:end-1)];
data_right = [data(:,2:end), zeros(rows,1)];
% Compute Laplacian using matrix operations
L = data_left + data_right - 2*data;
L_mag = abs(L);

% Auto-determine threshold if not provided
if isempty(threshold)
    threshold = prctile(L_mag(:), 95);
end

%% Step 2: Edge-Pair Detection
fprintf('Detecting edge pairs...\n');

% Initialize streak mask
streak_mask = false(size(data));

% For each row, find edge pairs
for i = 1:rows
    % Get Laplacian values for this row
    L_row = L(i, :);
    L_mag_row = L_mag(i, :);
    
    % Find all points above threshold
    edge_indices = find(L_mag_row >= threshold);
    
    % Look for edge pairs with opposite signs within MaxStreakWidth
    for j = 1:length(edge_indices)
        left_edge = edge_indices(j);
        
        % Look for right edge within MaxStreakWidth
        for k = (j+1):length(edge_indices)
            right_edge = edge_indices(k);
            
            % Check if edges are within MaxStreakWidth
            if (right_edge - left_edge) <= max_width
                % Check if edges have opposite signs
                if sign(L_row(left_edge)) ~= sign(L_row(right_edge))
                    % Mark the region between edges as streak
                    streak_mask(i, left_edge:right_edge) = true;
                end
            else
                % Edges are too far apart, skip remaining pairs
                break;
            end
        end
    end
end

%% Step 3: Connected Component Analysis
fprintf('Analyzing connected components...\n');

% Find connected components
CC = bwconncomp(streak_mask, 8); % 8-connectivity

% Initialize component information
streak_components = struct('PixelIdxList', {}, 'BoundingBox', {}, 'Width', {}, 'Height', {}, 'IsStreak', {});

for k = 1:CC.NumObjects
    pixels = CC.PixelIdxList{k};
    [rows_comp, cols_comp] = ind2sub(size(data), pixels);
    
    % Calculate component properties
    min_col = min(cols_comp);
    max_col = max(cols_comp);
    min_row = min(rows_comp);
    max_row = max(rows_comp);
    
    width = max_col - min_col + 1;
    height = max_row - min_row + 1;
    
    % Store component information
    component_info = struct();
    component_info.PixelIdxList = pixels;
    component_info.BoundingBox = [min_col, min_row, width, height];
    component_info.Width = width;
    component_info.Height = height;
    component_info.IsStreak = true; % All detected components are valid streaks
    
    streak_components(k) = component_info;
end

fprintf('Found %d streak components\n', length(streak_components));

%% Step 4: Visualization
visualization_data = struct();
visualization_data.original_data = data;
visualization_data.laplacian = L;
visualization_data.laplacian_magnitude = L_mag;
visualization_data.final_mask = streak_mask;
visualization_data.components = streak_components;
visualization_data.threshold = threshold;
visualization_data.max_width = max_width;

if show_viz
    createStreakVisualization(visualization_data);
end

end

function createStreakVisualization(viz_data)
    % Create comprehensive visualization
    figure('Name', 'Edge-Pair Streak Detection Analysis', 'Position', [100, 100, 1400, 1000]);
    
    % Original data
    subplot(2, 3, 1);
    imagesc(viz_data.original_data);
    title('Original Data');
    axis square;
    colormap parula;
    colorbar;
    
    % Laplacian
    subplot(2, 3, 2);
    imagesc(viz_data.laplacian);
    title('Laplacian');
    axis square;
    colormap parula;
    colorbar;
    
    % Laplacian magnitude with threshold
    subplot(2, 3, 3);
    imagesc(viz_data.laplacian_magnitude);
    title(sprintf('Laplacian Magnitude (Threshold: %.3f)', viz_data.threshold));
    axis square;
    colormap parula;
    colorbar;
    hold on;
    [edge_rows, edge_cols] = find(viz_data.laplacian_magnitude >= viz_data.threshold);
    plot(edge_cols, edge_rows, 'r.', 'MarkerSize', 1);
    hold off;
    
    % Final streak mask
    subplot(2, 3, 4);
    imagesc(viz_data.final_mask);
    title('Detected Streaks');
    axis square;
    colormap gray;
    colorbar;
    
    % Overlay on original data
    subplot(2, 3, 5);
    overlay_data = viz_data.original_data;
    overlay_data(viz_data.final_mask) = max(overlay_data(:)); % Highlight streaks
    imagesc(overlay_data);
    title('Streaks Overlaid on Original Data');
    axis square;
    colormap parula;
    colorbar;
    
    % Component analysis
    subplot(2, 3, 6);
    component_image = zeros(size(viz_data.original_data));
    
    for k = 1:length(viz_data.components)
        component_image(viz_data.components(k).PixelIdxList) = k;
    end
    
    imagesc(component_image);
    title('Streak Components (Colored)');
    axis square;
    colormap(lines);
    colorbar;
    
    % Add statistics
    valid_components = length(viz_data.components);
    streak_pixels = sum(viz_data.final_mask(:));
    coverage = 100 * streak_pixels / numel(viz_data.final_mask);
    
    fprintf('\nDetection Statistics:\n');
    fprintf('Threshold: %.3f\n', viz_data.threshold);
    fprintf('Max Streak Width: %d\n', viz_data.max_width);
    fprintf('Total Components: %d\n', valid_components);
    fprintf('Streak Pixels: %d\n', streak_pixels);
    fprintf('Coverage: %.1f%%\n', coverage);
    
    if valid_components > 0
        fprintf('\nStreak Component Details:\n');
        fprintf('Comp\tWidth\tHeight\n');
        for k = 1:length(viz_data.components)
            fprintf('%d\t%d\t%d\n', k, viz_data.components(k).Width, viz_data.components(k).Height);
        end
    end
end