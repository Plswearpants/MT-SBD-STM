function [similarity_score, qpi_a, qpi_b] = fftspace_corr2(kernel_a, kernel_b)
% Compute corr2 similarity in FFT-magnitude (Q-space).
%
% Inputs:
%   kernel_a, kernel_b : numeric 2D arrays
%
% Outputs:
%   similarity_score : corr2 between normalized Q-space magnitudes
%   qpi_a            : normalized abs(fftshift(fft2(kernel_a)))
%   qpi_b            : normalized abs(fftshift(fft2(kernel_b)))

    validateattributes(kernel_a, {'numeric'}, {'2d', 'nonempty'}, mfilename, 'kernel_a');
    validateattributes(kernel_b, {'numeric'}, {'2d', 'nonempty'}, mfilename, 'kernel_b');

    if ~isequal(size(kernel_a), size(kernel_b))
        error('fftspace_corr2:SizeMismatch', 'kernel_a and kernel_b must have identical size.');
    end

    qpi_a = abs(fftshift(fft2(kernel_a)));
    qpi_b = abs(fftshift(fft2(kernel_b)));

    max_a = max(qpi_a(:));
    max_b = max(qpi_b(:));
    if max_a <= 0 || max_b <= 0
        similarity_score = NaN;
        return;
    end

    qpi_a = qpi_a / max_a;
    qpi_b = qpi_b / max_b;
    similarity_score = corr2(qpi_a, qpi_b);
end
