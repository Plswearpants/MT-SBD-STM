% Create a 2D grid for x and y
[x, y] = meshgrid(-2:0.05:2, -2:0.05:2);

% --- 1. Non-convex energy landscape
E1 = sin(3*x).*cos(3*y) + 0.1*(x.^2 + y.^2);

% --- 2. Convexified energy landscape
lambda = 2.0;  % regularization strength
E2 = E1 + lambda * (x.^2 + y.^2);  % L2 convexification

% --- Plotting as 3D surfaces
figure;

subplot(1, 2, 1)
surf(x, y, E1, 'EdgeColor', 'none')  % smooth surface
title('Non-Convex Energy Landscape')
xlabel('x'), ylabel('y'), zlabel('Energy')
colormap turbo
view(45, 30)  % 3D view angle
colorbar
axis tight

subplot(1, 2, 2)
surf(x, y, E2, 'EdgeColor', 'none')
title('Convexified Energy Landscape')
xlabel('x'), ylabel('y'), zlabel('Energy')
colormap turbo
view(45, 30)
colorbar
axis tight
