%% Block 1: Load the .3ds data

% INPUTS
% 1: Data file to load, including file type ('QPI.3ds' for example)
% 2: Smoothing sigma for current data

% OUTPUTS
% header: Variable containing all experimental parameters
% I: Current data, smoothed by sigma
% dIdV: Numerically differentiated current data
% voltage: Vector of voltages for current
% midV: Vector on voltages for dIdV/QPI (midpoint of voltage vector)
% QPI: Fourier transformed dIdV data

% Modified function load3dsall from supplied matlab code from Nanonis
[header, par, I, dIdV, LockindIdV, bias, midV, QPI, LockinQPI] = load3dsall('Grid Spectroscopy008.3ds', 5);
xsize = header.grid_dim(1);
ysize = header.grid_dim(2);
elayer = header.points;
estart = par(1);
eend = par(2);
energy_range = linspace(estart, eend, elayer);
data_original = dIdV;
num_slices = size(data_original,3);
spatial = size(data_original,1);
rangetype = 'dynamic';

%% Original data view 
data_carried = data_original;
rangetype='dynamic';
figure;
d3gridDisplay(data_carried,rangetype);
slice_normalize = input('slice to normalize: ');
[data_carried] = normalizeBackgroundToZeroMean3D(data_carried, rangetype, slice_normalize);
%% remove bad columns
list_remove = [48,49];
data_carried(:,list_remove,:)= 0; % can put list_remove in the middle dimension.
Left = data_carried(:,1:47,:);
Right = data_carried(:,50:132,:);

%% variance level
figure; 
d3gridDisplay(data_carried, rangetype);
slice_num = input('which slice to use?');

[normalized_Left_3d, normalized_Right_3d, ~, ~] = equalVarianceLeftRight(Left, Right, rangetype, slice_num);
data_leveled =zeros(size(data_carried));
data_leveled(:,50:132,:)=normalized_Right_3d;
data_leveled(:,1:47,:)=normalized_Left_3d;
data_carried = data_leveled;

%% Local streak removal and interpolation 
figure; 
d3gridDisplay(data_carried, rangetype);
ref_slice = input('reference streaks slices: ');
[~, ~, streak_indices] = interpolateLocalStreaks(data_carried, ref_slice, []);
%% Prolifeartion streak removal for all slices at the same indices
data_corrected = zeros(size(data_carried));
for s = 1:size(data_carried,3)
    [data_corrected(:,:,s)] = interpolateLocalStreaks(data_carried, s, [], streak_indices);
end
data_carried = data_corrected;

%% Global streak correction
% Normalize background 
[data_carried] = normalizeBackgroundToZeroMean3D(data_carried, rangetype, ref_slice);

[data_streakremoved, QPI_nostreaks] = RemoveStreaks(data_carried, 'Direction', 'vertical');
data_carried = data_streakremoved;
%% Surface plot buffet
figure; 
surf(data_carried(:,:,125));

%%
function [normalized_Left_3d, normalized_Right_3d, background_info, comment] = equalVarianceLeftRight(Left, Right, rangeType, selected_slice)
% Normalizes left and right parts of a 3D dataset to have equal variance.
%   This function allows the user to interactively select background regions from both
%   left and right parts of a slice in a 3D dataset. The variance of the selected areas
%   is used to normalize the data so that both parts have equal variance.
%
% Arguments:
%   Left           3D array containing the left part of the data to be processed.
%   Right          3D array containing the right part of the data to be processed.
%   rangeType      Type of range for visualization ('global' or 'dynamic').
%   selected_slice (Optional) The index of the slice to use for background selection. Default is 1.
%
% Returns:
%   normalized_Left_3d   The normalized left part of the 3D data array.
%   normalized_Right_3d  The normalized right part of the 3D data array.
%   background_info      [left_variance, right_variance] of the selected background areas.
%   comment              Comment for logging the function call.
%
% Example:
%   [norm_Left, norm_Right, info, comment] = equalVarianceLeftRight(Left, Right, 'dynamic');
%   This example normalizes the left and right parts of the dataset based on the background
%   variance from the first slice.

arguments
    Left
    Right
    rangeType
    selected_slice = 1  % Default to slice 1 if not provided
end

% LOG comment of function call
comment = sprintf("equalVarianceLeftRight(selected_slice:%d)|", selected_slice);

% Get the dimensions of the input arrays
[rows_left, cols_left, k] = size(Left);
[rows_right, cols_right, ~] = size(Right);

% Display the selected slice and allow the user to select background regions
figure('Name', 'Select Background Areas');
subplot(1,2,1);
imagesc(Left(:,:,selected_slice));
title('Select Left Background Area');
axis square;
h_left = imrect;
position_left = wait(h_left);

subplot(1,2,2);
imagesc(Right(:,:,selected_slice));
title('Select Right Background Area');
axis square;
h_right = imrect;
position_right = wait(h_right);
close;

% Get the coordinates of the selected rectangles
x1_left = round(position_left(1));
y1_left = round(position_left(2));
x2_left = round(position_left(1) + position_left(3));
y2_left = round(position_left(2) + position_left(4));

x1_right = round(position_right(1));
y1_right = round(position_right(2));
x2_right = round(position_right(1) + position_right(3));
y2_right = round(position_right(2) + position_right(4));

% Initialize output arrays
normalized_Left_3d = zeros(size(Left));
normalized_Right_3d = zeros(size(Right));

% Process each slice
for i = 1:k
    % Extract background areas
    left_background = Left(y1_left:y2_left, x1_left:x2_left, i);
    right_background = Right(y1_right:y2_right, x1_right:x2_right, i);
    
    % Calculate variances
    left_variance = var(left_background(:));
    right_variance = var(right_background(:));
    
    % Store background info
    background_info = [left_variance, right_variance];
    
    % Calculate normalization ratio
    ratio = sqrt(right_variance / left_variance);
    
    % Normalize the data
    normalized_Left_3d(:,:,i) = Left(:,:,i) * ratio;
    normalized_Right_3d(:,:,i) = Right(:,:,i);
end

end
