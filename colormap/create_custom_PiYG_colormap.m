function create_custom_PiYG_colormap()
%CREATE_CUSTOM_PIYG_COLORMAP Creates and sets a custom colormap
%   This colormap is designed to emphasize increasing value,
%   with 3 distinct yet smoothly transitioning regions:
%     - [0, 0.85]: pure red to gold
%     - [0.85, 0.98]: gold to lime
%     - [0.98, 1.0]: lime to teal-green

    % Define normalized positions and RGB triplets
    positions = [0.00; 0.85; 0.98; 1.00];
    hex_colors = {'#FF0000', '#FFD700', '#DCE93D', '#00FF00'};

    % Convert hex colors to RGB (force row vector output)
    rgb_colors = zeros(length(hex_colors), 3);
    for i = 1:length(hex_colors)
        rgb_colors(i,:) = hex2rgb(hex_colors{i});
    end

    % Number of colors in the colormap
    n_colors = 256;
    x = linspace(0, 1, n_colors)'; % COLUMN vector

    % Interpolate RGB channels
    r = interp1(positions, rgb_colors(:,1), x, 'linear');
    g = interp1(positions, rgb_colors(:,2), x, 'linear');
    b = interp1(positions, rgb_colors(:,3), x, 'linear');

    custom_map = [r, g, b];  % For compatibility with previous code

    % Apply the colormap
    colormap(custom_map);
    colorbar;
    title('Custom PiYG-like Colormap');

    % Save the colormap for compatibility
    save('custom_PiYG_colormap.mat', 'custom_map');

end

function rgb = hex2rgb(hex)
%HEX2RGB Converts hex color string to RGB triplet in [0,1]
    hex = char(hex);
    if hex(1) == '#'
        hex = hex(2:end);
    end
    rgb = [hex2dec(hex(1:2)), hex2dec(hex(3:4)), hex2dec(hex(5:6))] / 255;
    rgb = reshape(rgb, 1, 3); % Ensure row vector
end 