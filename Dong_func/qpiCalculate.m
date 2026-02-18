function [QPI, comment] = qpiCalculate(data, m)
%This function calculates the Quasiparticle Interference (QPI) patterns using
%Fourier transform of input data (typically dI/dV or Lockin dI/dV)
%
% Dong Chen 2025/04
%
%   Input:
%   data        this is a 3d array form (x, y, V)
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
    validateattributes(m, {'numeric'}, {'scalar','integer','positive','finite','real'});
    QPI = zeros(m, m, points);
    calculateQPI = @(slice) abs(fftshift(fft2(slice - mean(mean(slice)), m, m)));
end

for i = 1:points
    QPI(:,:,i) = calculateQPI(data(:,:,i));
end

end
