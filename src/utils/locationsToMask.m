function mask = locationsToMask(locations, mask_size)
% LOCATIONSTOMASK Convert a list of locations to a binary mask
%   mask = locationsToMask(locations, mask_size) creates a binary mask where
%   points closest to the given locations are 1 and others are 0
%
% Inputs:
%   locations - Nx2 array of [x,y] coordinates (can be non-integer)
%   mask_size - [height, width] of the desired mask
%
% Outputs:
%   mask - Binary mask of size mask_size
%
% Example:
%   locations = [10.5, 20.3; 30.2, 40.1];
%   mask = locationsToMask(locations, [100, 100]);

% Create coordinate grids
[X, Y] = meshgrid(1:mask_size(2), 1:mask_size(1));

mask = zeros(mask_size);
% Calculate minimum distance to any location for each point
for i = 1:size(locations, 1)
    % Calculate distance from each point to current location
    dist = sqrt((X - locations(i,1)).^2 + (Y - locations(i,2)).^2);
    % Update minimum distances
    min_d=min(dist,[],'all');
    [x,y]=find(dist==min_d);
    mask(x,y)=1;
end
end 