function [D, D_magnitude] = computeSpatialDerivative(matrix, direction)
% COMPUTESPATIALDERIVATIVE Computes spatial derivative of a 2D matrix in specified direction
%   Inputs:
%       matrix: 2D input matrix
%       direction: 'x' or 'y' for horizontal or vertical derivative
%   Outputs:
%       D: Difference map in specified direction
%       D_magnitude: Magnitude of the difference map
%
%   Example:
%       [D, D_magnitude] = computeSpatialDerivative(myMatrix, 'x');

% Input validation
if ~ismatrix(matrix)
    error('Input must be a 2D matrix');
end

if ~ismember(direction, {'x', 'y'})
    error('Direction must be either ''x'' or ''y''');
end

% Compute difference map
if strcmp(direction, 'x')
    % Horizontal derivative (differences along columns)
    D = diff(matrix, 1, 2);
    % Pad with zeros to maintain size
    D = [D, zeros(size(matrix, 1), 1)];
else
    % Vertical derivative (differences along rows)
    D = diff(matrix, 1, 1);
    % Pad with zeros to maintain size
    D = [D; zeros(1, size(matrix, 2))];
end

% Compute magnitude of difference map
D_magnitude = abs(D);

% Visualize results
figure('Name', 'Spatial Derivative Analysis');

% Original matrix
subplot(1, 3, 1);
imagesc(matrix);
title('Original Matrix');
colorbar;
axis square;

% Difference map
subplot(1, 3, 2);
imagesc(D);
title(['Difference Map (', direction, '-direction)']);
colorbar;
axis square;

% Magnitude of difference map
subplot(1, 3, 3);
imagesc(D_magnitude);
title(['Magnitude of Difference Map (', direction, '-direction)']);
colorbar;
axis square;

% Adjust colormap for better visualization
colormap('jet');
end 