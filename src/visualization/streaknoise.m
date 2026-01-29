Y_in = data_carried(:,:,65);
figure; imagesc(Y_in); axis square; colormap('gray')
figure; surf(Y_in(40:50,:)); axis square;
%%
% Create a 128x128 dataset with sinusoidal pattern
[X,Y] = meshgrid(1:128, 1:128);
freq = 0.3; % Frequency of the sine wave
amplitude = 1; % Amplitude of the sine wave
Z = amplitude * sin(freq * (X+Y));

figure; imagesc(Z); axis square; colormap('gray'); title('Sinusoidal Pattern');
figure; surf(Z, 'EdgeColor', 'none'); axis square; title('3D View of Sinusoidal Pattern');

%%
Z(:,45) = Z(:,45)+1;
Z(:,70:72) = Z(:,70:72)-0.5;
figure; surf(Z, 'EdgeColor', 'none'); axis square; title('3D View of Sinusoidal Pattern');
figure; imagesc(Z); axis square; colormap('gray'); title('Sinusoidal Pattern');
