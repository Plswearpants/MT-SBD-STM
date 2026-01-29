function [adjusted_data, noise_increase_std, noise_std, comment] = levelNoiseInteractive(data, dimension)
% Levels noise across a specified dimension in a 2D dataset by increasing the noise variance.
%   This function allows the user to interactively select a region containing 
%   noise in a 2D dataset, then levels the noise across the specified dimension.
%   The noise variance is increased in lower noise regions to match the highest noise variance,
%   and the noise map is adjusted to maintain the standard deviation relative to the mean.
%
% Arguments:
%   data        2D array containing the data to be processed.
%   dimension   Dimension along which to level the noise ('x' for rows, 'y' for columns).
%
% Returns:
%   adjusted_data    The data array with noise leveled across the specified dimension.
%   noise_increase   The amount of noise variance added along the selected dimension.
%   comment          Comment for logging the function call.
%
% August 2024 - Dong Chen
%
% Example:
%   [adjusted_data, noise_increase, comment] = levelNoiseInteractive(data, 'x');
%   This example levels the noise variance along the x-dimension interactively.

arguments
    data
    dimension
end

% LOG comment of function call
comment = sprintf("levelNoiseInteractive(dimension:%s)|", dimension);

% Display the data and let the user draw the window W
figure;
imagesc(data);  % Display the data as an image
colormap('gray');
title('Select the window containing the noise to level');
h = drawrectangle('Color','r');  % User draws the rectangle
position = round(h.Position);    % Get the position of the rectangle [x, y, width, height]

% Process the selected window based on the specified dimension
if dimension == 'x'
    % Extend W along rows
    W_start = [position(2), position(1)];
    W_end = [position(2) + position(4) - 1, position(1) + position(3) - 1];
    W_mod = data(W_start(1):W_end(1), :);  % Extend W to the full width of the data
    noise_std = std(W_mod, 0, 1);          % Calculate variance along columns in W_mod

elseif dimension == 'y'
    % Extend W along columns
    W_start = [position(2), position(1)];
    W_end = [position(2) + position(4) - 1, position(1) + position(3) - 1];
    W_mod = data(:, W_start(2):W_end(2));  % Extend W to the full height of the data
    noise_std = std(W_mod, 0, 2);          % Calculate variance along rows in W_mod

else
    error('Invalid dimension. Choose either "x" or "y".');
end

% Find the maximum noise variance and mean in W_mod
max_noise_std = max(noise_std);

% Calculate the required increase in noise variance
noise_increase_std = max_noise_std - noise_std;

% Generate a random matrix with entries between 0 and 1
if dimension == 'x'
    noise_map = randn(size(data, 1), size(data, 2)); % Random zero-mean gaussian noise
    noise_map = noise_map .* (repmat(noise_increase_std, size(data, 1), 1)); 
elseif dimension == 'y'
    noise_map = randn(size(data, 1), size(data, 2)); % Random zero-mean gaussian noise
    noise_map = noise_map .* (repmat(noise_increase_std', 1, size(data, 2)));
end

% Adjust the data by adding the generated noise map
adjusted_data = data + noise_map;

% Optional: Plot the adjusted data (can be commented out if not needed)
figure;
imagesc(adjusted_data);
colormap('gray');
colorbar;
title('Adjusted Data with Leveled Noise');
axis image;

end
