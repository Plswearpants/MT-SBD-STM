function [QPI, comment] = qpiCalculate(data, m)
%This function calculates the Quasiparticle Interference (QPI) patterns using
%Fourier transform of input data (typically dI/dV or Lockin dI/dV)
%
% Dong Chen 2025/04
%
%   Input:
%   data        this is a 3d array form (x, y, V)
%   m           FFT output size:
%               - scalar m -> use square size (m, m)
%               - vector [m, n] -> use non-square size (m, n)
%
%   Output:
%   QPI         this is a 3d array form (x, y, V) in the fourier domain 

arguments
    data
    m = []
end

[x_num, y_num, points] = size(data);

comment = sprintf("qpiCalculate(data:%s x %s x %s)|", ...
    mat2str(x_num), mat2str(y_num), mat2str(points));

if isempty(m)
    QPI = zeros(x_num, y_num, points);
    calculateQPI = @(slice) abs(fftshift(fft2(slice - mean(mean(slice)))));
else
    validateattributes(m, {'numeric'}, {'vector','integer','positive','finite','real'});
    if numel(m) == 1
        m_out = m;
        n_out = m;
    elseif numel(m) == 2
        m_out = m(1);
        n_out = m(2);
    else
        error('qpiCalculate:InvalidSizeInput', ...
            'm must be empty, scalar, or a 2-element vector [m, n].');
    end
    QPI = zeros(m_out, n_out, points);
    calculateQPI = @(slice) abs(fftshift(fft2(slice - mean(mean(slice)), m_out, n_out)));
end

for i = 1:points
    QPI(:,:,i) = calculateQPI(data(:,:,i));
end

end
