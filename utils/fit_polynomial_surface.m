function [fitted_surface, subtracted_surface] = fit_polynomial_surface(data, degree, varargin)
% FIT_POLYNOMIAL_SURFACE Fits a polynomial surface to data with optional polygon masks
%   [fitted_surface, subtracted_surface] = fit_polynomial_surface(data, degree, ...)
%
% Inputs:
%   data: 2D or 3D array of data points
%   degree: [x_deg, y_deg] array specifying polynomial degrees for x and y directions
%   Optional parameters:
%       'Mask': Polygon mask(s) for fitting (cell array of Nx2 points)
%               If not provided, interactive mask drawing will be enabled
%       'Order': Order of polynomial terms (default: 'full')
%       'Vis': Visualization flag (default: 0)
%              If 1, shows original data, fitted surface, and subtracted surface
%
% Outputs:
%   fitted_surface: The fitted polynomial surface
%   subtracted_surface: Original data minus fitted surface
%
% Example:
%   data = peaks(100);
%   % Interactive mask drawing
%   [fit, sub] = fit_polynomial_surface(data, [2 3], 'Vis', 1);
%   % Or with predefined masks
%   mask = {[20 20; 20 80; 80 80; 80 20]}; % Square mask
%   [fit, sub] = fit_polynomial_surface(data, [2 3], 'Mask', mask, 'Vis', 1);

% Parse input arguments
p = inputParser;
addRequired(p, 'data', @(x) isnumeric(x) && (ndims(x) == 2 || ndims(x) == 3));
addRequired(p, 'degree', @(x) isnumeric(x) && isvector(x) && length(x) == 2 && all(x >= 0));
addParameter(p, 'Mask', {}, @(x) iscell(x) || isempty(x));
addParameter(p, 'Order', 'full', @(x) ischar(x) && ismember(x, {'full', 'reduced'}));
addParameter(p, 'Vis', 0, @(x) isnumeric(x) && isscalar(x) && ismember(x, [0 1]));
parse(p, data, degree, varargin{:});

% Get parameters
data = p.Results.data;
degree_x = degree(1);
degree_y = degree(2);
masks = p.Results.Mask;
order = p.Results.Order;
vis_flag = p.Results.Vis;

% Handle interactive mask drawing if no masks provided
if isempty(masks)
    masks = draw_interactive_masks(data);
end

% Handle 3D data
if ndims(data) == 3
    [ny, nx, nz] = size(data);
    fitted_surface = zeros(size(data));
    subtracted_surface = zeros(size(data));
    
    % Process each slice
    for z = 1:nz
        [fit_slice, sub_slice] = fit_single_surface(data(:,:,z), degree_x, degree_y, masks, order);
        fitted_surface(:,:,z) = fit_slice;
        subtracted_surface(:,:,z) = sub_slice;
    end
    
    % Visualization for 3D data (show middle slice)
    if vis_flag
        z_mid = round(nz/2);
        visualize_results(data(:,:,z_mid), fitted_surface(:,:,z_mid), subtracted_surface(:,:,z_mid), masks);
    end
    return;
end

% For 2D data, process single surface
[fitted_surface, subtracted_surface] = fit_single_surface(data, degree_x, degree_y, masks, order);

% Visualization for 2D data
if vis_flag
    visualize_results(data, fitted_surface, subtracted_surface, masks);
end

end

function masks = draw_interactive_masks(data)
    % Create figure for mask drawing
    figure('Name', 'Draw Masks (Press Enter when done)');
    imagesc(data);
    colormap(gray);
    axis equal tight;
    title('Draw masks (Press Enter when done)');
    
    masks = {};
    while true
        % Draw polygon
        h = drawpolygon('Color', 'r');
        
        % Check if user pressed Enter
        if isempty(h.Position)
            break;
        end
        
        % Add mask to cell array
        masks{end+1} = h.Position;
        
        % Change color for next mask
        h.Color = 'g';
    end
end

function visualize_results(data, fitted, subtracted, masks)
    figure('Name', 'Surface Fitting Results', 'Position', [100 100 1200 400]);
    
    % Original data with masks
    subplot(1,3,1);
    imagesc(data);
    colormap(gray);
    axis equal tight;
    title('Original Data with Masks');
    hold on;
    for i = 1:length(masks)
        plot(masks{i}(:,1), masks{i}(:,2), 'r-', 'LineWidth', 2);
    end
    hold off;
    
    % Fitted surface
    subplot(1,3,2);
    imagesc(fitted);
    colormap(gray);
    axis equal tight;
    title('Fitted Surface');
    
    % Subtracted surface
    subplot(1,3,3);
    imagesc(subtracted);
    colormap(gray);
    axis equal tight;
    title('Subtracted Surface');
end

function [fitted_surface, subtracted_surface] = fit_single_surface(data, degree_x, degree_y, masks, order)
    [ny, nx] = size(data);
    [X, Y] = meshgrid(1:nx, 1:ny);
    
    % Create mask if none provided
    if isempty(masks)
        mask = true(ny, nx);
    else
        mask = false(ny, nx);
        for i = 1:length(masks)
            poly_mask = poly2mask(masks{i}(:,1), masks{i}(:,2), ny, nx);
            mask = mask | poly_mask;
        end
    end
    
    % Get points within mask
    x = X(mask);
    y = Y(mask);
    z = data(mask);
    
    % Generate polynomial terms based on order
    if strcmp(order, 'full')
        terms = generate_polynomial_terms(degree_x, degree_y);
    else
        terms = generate_reduced_terms(degree_x, degree_y);
    end
    
    % Create design matrix
    A = zeros(length(x), size(terms, 1));
    for i = 1:size(terms, 1)
        A(:,i) = x.^terms(i,1) .* y.^terms(i,2);
    end
    
    % Solve least squares problem
    coeffs = A \ z;
    
    % Generate fitted surface
    fitted_surface = zeros(ny, nx);
    for i = 1:size(terms, 1)
        fitted_surface = fitted_surface + coeffs(i) * X.^terms(i,1) .* Y.^terms(i,2);
    end
    
    % Calculate subtracted surface
    subtracted_surface = data - fitted_surface;
end

function terms = generate_polynomial_terms(degree_x, degree_y)
    terms = [];
    for i = 0:degree_x
        for j = 0:degree_y
            terms = [terms; i j];
        end
    end
end

function terms = generate_reduced_terms(degree_x, degree_y)
    terms = [];
    for i = 0:degree_x
        for j = 0:degree_y
            if i == 0 || j == 0 || i == j
                terms = [terms; i j];
            end
        end
    end
end 