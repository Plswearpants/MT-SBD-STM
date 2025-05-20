function cmap = invgray(m)
%INVGRAY Inverted grayscale colormap
%   INVGRAY(M) returns an M-by-3 matrix containing an inverted grayscale colormap.
%   INVGRAY, by itself, is the same length as the current figure's colormap.
%   If no figure exists, MATLAB uses the length of the default colormap.
%
%   Example:
%       colormap(invgray)
%       colormap(invgray(128))
%
%   See also GRAY, COLORMAP, RGBPLOT.

if nargin < 1
   f = get(groot,'CurrentFigure');
   if isempty(f)
      m = size(get(groot,'DefaultFigureColormap'),1);
   else
      m = size(f.Colormap,1);
   end
end

% Create inverted grayscale colormap
cmap = flipud(gray(m));
end 