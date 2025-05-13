function [data_cleaned, fit_params] = lorentzianBraggRemove(data, varargin)
% LORENTZIANBRAGGREMOVE Removes Bragg peaks using asymmetric Lorentzian fitting
%   This function removes Bragg peaks from 2D data using asymmetric Lorentzian
%   fitting in both x and y directions with Gaussian broadening.
%   The peak selection is done in q-space (Fourier space) while input and output
%   are in real space.
%
% Input:
%   data: 2D array containing the data with Bragg peaks (real space)

%
% Output:
%   data_cleaned: data with Bragg peaks removed (real space)
%   fit_params: structure containing fitting parameters for each peak
%
% Example:
%   [data_cleaned, fit_params] = lorentzianBraggRemove(data);
%   [data_cleaned, fit_params] = lorentzianBraggRemove(data, 'sigma', 1.5);
%
% See also: asymmetric_lorentzian, fit_lorentzian

% Initialize
data_cleaned = data;
fit_params = struct();
iteration = 1;

% Convert to q-space for peak selection
QPI = fftshift(fft2(data_cleaned));
QPI_logabs = log(abs(QPI));

% Create figure for interactive peak selection in q-space
fig = figure('Name', 'Select Bragg Peak Region in Q-space');
imagesc(QPI_logabs);
colormap(gray);
colorbar;
axis square;
title('Select Bragg peaks in Q-space');

while true
    % 1. Let user draw a square around the Bragg peak in q-space
    title(sprintf('Draw square around Bragg peak %d in Q-space', iteration));
    h = drawrectangle('Color', 'r');
    wait(h);
    
    % Get rectangle parameters
    rect = h.Position;  % [x y width height]
    x1 = round(rect(1));
    y1 = round(rect(2));
    x2 = round(x1 + rect(3));
    y2 = round(y1 + rect(4));
    
    % 2. Find highest intensity point in the square in q-space
    [x, y] = meshgrid(1:size(QPI,2), 1:size(QPI,1));
    mask = (x >= x1 & x <= x2 & y >= y1 & y <= y2);
    [max_val, max_idx] = max(abs(QPI(mask)));
    [peak_y, peak_x] = ind2sub(size(QPI), find(mask));
    peak_x = peak_x(max_idx);
    peak_y = peak_y(max_idx);
    
    % 3. Extract data within the square for fitting
    x_data = abs(QPI(peak_y, x1:x2));  % Horizontal line through peak
    y_data = abs(QPI(y1:y2, peak_x));  % Vertical line through peak
    
    % Store original data for later use
    x_data_orig = x_data;
    y_data_orig = y_data;
    
    % Normalize the data
    x_data_norm = x_data / max(x_data);
    y_data_norm = y_data / max(y_data);
    
    % Create position vectors for fitting
    x_pos = x1:x2;
    y_pos = y1:y2;
    
    % Fit x-direction (horizontal line)
    x_fit = fit_lorentzian(x_pos, x_data_norm, peak_x);
    
    % Fit y-direction (vertical line)
    y_fit = fit_lorentzian(y_pos, y_data_norm, peak_y);
    
    % 4. Create 2D model and fit sigma
    [X, Y] = meshgrid(1:size(QPI,2), 1:size(QPI,1));
    
    % Get the normalization factors
    x_norm_factor = max(x_data_orig);
    y_norm_factor = max(y_data_orig);
    
    % Create 1D profiles
    x_profile = x_fit(X(1,:));
    y_profile = y_fit(Y(:,1));
    
    % Create initial 2D model
    model_norm = y_profile * x_profile;
    model = model_norm * x_norm_factor * y_norm_factor;
    
    % Extract the region around the peak for sigma fitting
    peak_region = abs(QPI(y1:y2, x1:x2));
    
    % Define the objective function for sigma fitting
    objective_func = @(sigma) fit_sigma(sigma, model(y1:y2, x1:x2), peak_region);
    
    % Fit sigma using fminbnd
    sigma = fminbnd(objective_func, 0.1, 5);
    
    % Apply the fitted sigma
    [Xg, Yg] = meshgrid(-3*sigma:3*sigma, -3*sigma:3*sigma);
    gaussian_kernel = exp(-(Xg.^2 + Yg.^2)/(2*sigma^2));
    gaussian_kernel = gaussian_kernel / sum(gaussian_kernel(:));
    
    % Apply 2D Gaussian broadening
    model = conv2(model, gaussian_kernel, 'same');
    
    % Scale the model to match the max of the data in the selected region
    model_region = model(y1:y2, x1:x2);
    data_region = abs(QPI(y1:y2, x1:x2));
    scale_factor = max(data_region(:)) / max(model_region(:));
    model = model * scale_factor;

    % 5. Remove the peak in q-space
    QPI_cleaned = QPI - model;  % Subtractive removal, now with scaled model
    
    % Plot the original, model, and cleaned square regions only
    fit_fig = figure('Name', sprintf('Bragg Peak Removal Visualization - Peak %d', iteration));
    
    % Original square region
    subplot(1,3,1);
    orig_data = (abs(QPI(y1:y2, x1:x2)));
    imagesc(abs(orig_data));
    colormap(gray);
    colorbar;
    axis square;
    title('Original Square Region');
    
    % Fitted model in the square region
    subplot(1,3,2);
    model_region = model(y1:y2, x1:x2);
    imagesc(model_region);
    colormap(gray);
    colorbar;
    axis square;
    title('Fitted 2D Model in Square');
    clim auto

    % Cleaned square region
    subplot(1,3,3);
    cleaned_data = abs(QPI_cleaned(y1:y2, x1:x2));
    imagesc(cleaned_data);
    colormap(gray);
    colorbar;
    axis square;
    title('Cleaned Square Region');
    
    % Use the same color scale for all images for direct comparison
    cmin = min([orig_data(:); model_region(:); cleaned_data(:)]);
    cmax = max([orig_data(:); model_region(:); cleaned_data(:)]);
    subplot(1,3,1); caxis([cmin cmax]);
    subplot(1,3,2); caxis([cmin cmax]);
    subplot(1,3,3); caxis([cmin cmax]);
    
    set(fit_fig, 'Position', [100, 100, 1500, 500]);
    
    % Convert back to real space
    data_cleaned = real(ifft2(ifftshift(QPI_cleaned)));
    
    % Store fit parameters and model
    fit_params(iteration).peak_position = [peak_x, peak_y];
    fit_params(iteration).x_params = x_fit;
    fit_params(iteration).y_params = y_fit;
    fit_params(iteration).sigma = sigma;
    fit_params(iteration).selection_rect = rect;
    fit_params(iteration).x_norm_factor = x_norm_factor;
    fit_params(iteration).y_norm_factor = y_norm_factor;
    fit_params(iteration).model = model;
    
    % Update the q-space display
    figure(fig);
    imagesc(log(abs(QPI_cleaned)));
    colormap(gray);
    colorbar;
    axis square;
    title(sprintf('Peak %d removed in Q-space', iteration));
    
    % Ask user if they want to continue
    choice = questdlg('Do you want to remove another peak?', ...
        'Continue?', ...
        'Yes', 'No', 'Yes');
    
    if strcmp(choice, 'No')
        break;
    end
    
    iteration = iteration + 1;
end

% Clean up
if ishandle(fig)
    close(fig);
end
end

function y = asymmetric_lorentzian(x, x0, gamma, alpha)
% ASYMMETRIC_LORENTZIAN Calculates asymmetric Lorentzian function
%   x0: peak position
%   gamma: width parameter
%   alpha: asymmetry parameter
    y = (1/pi) * (gamma/2) ./ ((x - x0).^2 + (gamma/2)^2) .* (1 + alpha*(x - x0));
end

function fit = fit_lorentzian(x, y, peak_pos)
% FIT_LORENTZIAN Fits asymmetric Lorentzian to 1D data
%   Returns function handle for the fit
    
    % Ensure x and y are column vectors
    x = x(:);
    y = y(:);
    
    % Initial guess for parameters
    x0 = peak_pos;
    gamma = 2;  % Initial width guess
    alpha = 0;  % Initial asymmetry guess
    
    % Define the model function
    model = @(params, x) asymmetric_lorentzian(x, params(1), params(2), params(3));
    
    % Define bounds for parameters
    lb = [peak_pos - 5, 0.1, -1];  % Lower bounds
    ub = [peak_pos + 5, 10, 1];    % Upper bounds
    
    % Fit using lsqcurvefit with bounds
    options = optimoptions('lsqcurvefit', ...
        'Display', 'off', ...
        'MaxFunctionEvaluations', 1000, ...
        'MaxIterations', 1000, ...
        'FunctionTolerance', 1e-6, ...
        'StepTolerance', 1e-6);
    
    params = lsqcurvefit(model, [x0, gamma, alpha], x, y, lb, ub, options);
    
    % Return function handle with fitted parameters
    fit = @(x) asymmetric_lorentzian(x, params(1), params(2), params(3));
end

function error = fit_sigma(sigma, model, data)
    % Create Gaussian kernel
    [Xg, Yg] = meshgrid(-3*sigma:3*sigma, -3*sigma:3*sigma);
    gaussian_kernel = exp(-(Xg.^2 + Yg.^2)/(2*sigma^2));
    gaussian_kernel = gaussian_kernel / sum(gaussian_kernel(:));
    
    % Apply Gaussian broadening
    model_broadened = conv2(model, gaussian_kernel, 'same');
    
    % Calculate error (using normalized data)
    data_norm = data / max(data(:));
    model_norm = model_broadened / max(model_broadened(:));
    error = sum((data_norm(:) - model_norm(:)).^2);
end