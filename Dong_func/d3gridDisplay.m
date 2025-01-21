function d3gridDisplay(LDoS_noisy, rangeType)
% Displays a 3D dataset using a specified colormap and range type.
%   This function visualizes a 3D dataset by converting each slice to an image 
%   using a colormap and either a global or dynamic range. The processed slices 
%   are then displayed as a 3D image stack.
%
% Arguments:
%   LDoS_noisy  3D array containing the data to be displayed.
%   rangeType   Type of range for visualization ('global' or 'dynamic').
%
% Returns:
%   None. The function displays the 3D dataset as an image stack.
%
% August 2024 - Dong Chen
%
% Example:
%   d3gridDisplay(LDoS_noisy, 'global');
%   This example displays the dataset with a global intensity range across slices.

% Load colormap
load('InverseGray', 'invgray');
map = invgray;

% Determine global min and max values if global range is selected
if strcmp(rangeType, 'global')
    globalMin = min(LDoS_noisy(:));
    globalMax = max(LDoS_noisy(:));
end

% Preallocate the sliced LDoS array
slicedglob_LDoS = zeros(size(LDoS_noisy, 1), size(LDoS_noisy, 2), size(LDoS_noisy, 3), 3);

% nos for dynamic range, larger nos bigger contrast
nos = 8;

% Convert LDoS_result to image format using the colormap
for k = 1:size(LDoS_noisy, 3)
    if strcmp(rangeType, 'dynamic')
        MeddIdV(k) = median(median(LDoS_noisy(:,:,k)));
        Stdv(k) = std(std(LDoS_noisy(:,:,k)));
        range = [MeddIdV(k)-nos*Stdv(k) MeddIdV(k)+nos*Stdv(k)];
    else
        range = [globalMin globalMax];
    end
    slicedglob_LDoS(:,:,k,:) = mat2im(LDoS_noisy(:,:,k), map, range);
end

% Display the 3D image stack
imshow3D(slicedglob_LDoS);

end
