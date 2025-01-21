function [B, selectedAreaMean, selectedAreaVariance, comment] = padWithSelectedAreaNoise(A, m, selected_slice)
% Pads each slice of a 3D dataset to a larger size with noise matching the statistics 
% of a user-selected area in the slice.
%   This function allows the user to interactively select an area in a slice of a 3D dataset.
%   The mean and variance of the selected area are used to generate Gaussian noise, 
%   which is then used to pad each slice of the dataset to the specified size.
%
% Arguments:
%   A               3D array containing the data to be processed (n*n*k).
%   m               Target dimension for the padded output (must be larger than n).
%   selected_slice  (Optional) The index of the slice to use for area selection. Default is 1.
%
% Returns:
%   B                       The padded 3D data array.
%   selectedAreaMean        The mean of the selected area.
%   selectedAreaVariance    The variance of the selected area.
%   comment                 Comment for logging the function call.
%
% August 2024 - Dong Chen
%
% Example:
%   [B, mean, variance, comment] = padWithSelectedAreaNoise(A, 10);
%   This example pads each slice of the dataset to 10x10 using the mean and variance
%   of the selected area in the first slice.

arguments
    A
    m
    selected_slice = 1  % Default to slice 1 if not provided
end

% LOG comment of function call
comment = sprintf("padWithSelectedAreaNoise(selected_slice:%d, target_size:%d)|", selected_slice, m);

% Get the dimensions of the input array A
[n, ~, k] = size(A);

% Initialize the output array B with zeros
B = zeros(m, m, k);

% Display the selected slice and get user input for the rectangle
slice = A(:,:,selected_slice);
figure, imshow(slice, []);
title(sprintf('Select an area for slice %d', selected_slice));
h = imrect;
position = wait(h);
close;

% Calculate the coordinates of the selected rectangle
x1 = round(position(1));
y1 = round(position(2));
x2 = round(position(1) + position(3));
y2 = round(position(2) + position(4));

% Extract the selected area
selectedArea = slice(y1:y2, x1:x2);

% Calculate the mean and variance of the selected area
selectedAreaMean = mean(selectedArea(:));
selectedAreaVariance = var(selectedArea(:));

% Loop over each slice in the third dimension
for i = 1:k
    % Extract the current slice
    currentSlice = A(:,:,i);
    
    % Create noise with the calculated mean and variance
    noise = selectedAreaMean + sqrt(selectedAreaVariance) * randn(m, m);
    
    % Calculate padding sizes
    padSizeX = [(m - n) / 2, (m - n) / 2];
    padSizeY = [(m - n) / 2, (m - n) / 2];
    
    if mod(m-n, 2) ~= 0
        padSizeX(2) = ceil(padSizeX(2));
        padSizeY(2) = ceil(padSizeY(2));
    end
    
    % Initialize the output slice with the noise
    paddedSlice = noise;
    
    % Place the original slice in the center
    paddedSlice(padSizeX(1)+1:padSizeX(1)+n, padSizeY(1)+1:padSizeY(1)+n) = currentSlice;
    
    % Replace the selected area in the padded slice with the original data
    startX = padSizeX(1) + x1;
    startY = padSizeY(1) + y1;
    endX = startX + (x2 - x1);
    endY = startY + (y2 - y1);
    paddedSlice(startY:endY, startX:endX) = currentSlice(y1:y2, x1:x2);
    
    % Store the result in the output array
    B(:,:,i) = paddedSlice;
end

end
