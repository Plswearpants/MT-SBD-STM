function create_custom_PiYG_colormap()
    % Create PiYG colormap with shifted divergent point
    n_colors = 256;
    
    % Define key colors from matplotlib PiYG
    colors = [
        142, 1, 82;     % Dark pink
        197, 27, 125;   % Pink
        222, 119, 174;  % Light pink
        241, 182, 218;  % Very light pink
        253, 224, 239;  % Near white pink
        247, 247, 247;  % White
        230, 245, 208;  % Near white green
        184, 225, 134;  % Very light green
        127, 188, 65;   % Light green
        77, 146, 33;    % Green
        39, 100, 25     % Dark green
    ] / 255;  % Normalize to [0,1]
    
    % Create interpolation points with shifted divergent point at 0.85
    n_points = size(colors, 1);
    mid_point = 6;  % Index of the white color (divergent point)
    
    % Create two sets of points: before and after divergent point
    x_before = linspace(0, 0.85, mid_point);  % Points up to 0.85
    x_after = linspace(0.85, 1, n_points - mid_point + 1);  % Points from 0.85 to 1
    x_points = [x_before(1:end-1) x_after];  % Combine, avoiding duplicate at 0.85
    
    % Create evenly spaced points for final colormap
    xi = linspace(0, 1, n_colors);
    
    % Interpolate colors
    r = interp1(x_points, colors(:,1), xi, 'pchip');
    g = interp1(x_points, colors(:,2), xi, 'pchip');
    b = interp1(x_points, colors(:,3), xi, 'pchip');
    
    custom_map = [r' g' b'];
    
    % Save the colormap
    save('custom_PiYG_colormap.mat', 'custom_map');
end 