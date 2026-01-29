% Create an inverted grayscale colormap
% This colormap goes from white (1) to black (0)

% Create a linear grayscale colormap
gray_map = linspace(1, 0, 256)';

% Create the RGB version (all channels are the same for grayscale)
invgray = repmat(gray_map, 1, 3);

% Save the colormap
save('invgray.mat', 'invgray');

% Display the colormap
figure;
colormap(invgray);
colorbar;
title('Inverted Grayscale Colormap'); 