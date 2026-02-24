function out = inplaneShift(img, anchor_from, anchor_to)
% Shift image so the point at anchor_from moves to anchor_to.
% Output has the same size as img. Pixels with no source are set to zero.
%
%   img          2D array (or 3D; shift applied per slice)
%   anchor_from  [row, col] of point to move (1-based)
%   anchor_to    [row, col] destination (1-based)

arguments
    img
    anchor_from (1,2) {mustBeNumeric,mustBeInteger,mustBePositive}
    anchor_to   (1,2) {mustBeNumeric,mustBeInteger,mustBePositive}
end

dr = anchor_to(1) - anchor_from(1);
dc = anchor_to(2) - anchor_from(2);

if ndims(img) == 2
    out = shiftOnce(img, dr, dc);
else
    [nr, nc, nz] = size(img);
    out = zeros(nr, nc, nz, 'like', img);
    for k = 1:nz
        out(:,:,k) = shiftOnce(img(:,:,k), dr, dc);
    end
end

end

function out = shiftOnce(img, dr, dc)
[nr, nc] = size(img);
[R, C] = ndgrid(1:nr, 1:nc);
src_r = R - dr;
src_c = C - dc;
valid = src_r >= 1 & src_r <= nr & src_c >= 1 & src_c <= nc;
out = zeros(nr, nc, 'like', img);
out(valid) = img(src_r(valid) + nr * (src_c(valid) - 1));
end
