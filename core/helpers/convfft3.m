function [ Y ] = convfft3( A, X, varargin )
%CONVFFT3   FFT implementation of 2D linear convolution for 3D arrays
%   Performs 2D convolution along the third dimension of 3D arrays
%
%   Y = convfft3(A, X) convolves each slice of A with X.
%   A must be a 3D array of kernels. X can be:
%   - 3D with matching third dimension (slice-wise convolution), or
%   - 2D / single-slice, which is applied to all slices of A.
%
%   Y = convfft3(..., adj, tmpsz, outsz)  
%   applies the adjoint operation if ADJ==TRUE
%   TMPSZ and OUTSZ sets the extension and output sizes manually for each slice
%
%   Example:
%       A = randn(10,10,3);  % 3D array of kernels
%       X = randn(20,20,3);  % 3D array of activations
%       Y = convfft3(A, X);  % Convolve each slice

    % Validate input dimensions
    if ndims(A) ~= 3
        error('A must be a 3D array.');
    end

    if ndims(X) > 3
        error('X must be a 2D or 3D array.');
    end

    if ismatrix(X) || size(X,3) == 1
        X = X(:,:,ones(1,size(A,3)));
    elseif size(A,3) ~= size(X,3)
        error('Third dimension of A and X must match, or X must be 2D/single-slice.');
    end
    
    % Parse optional arguments
    if numel(varargin) > 3
        error('Too many input arguments.');
    end
    
    if numel(varargin) >= 1 && ~isempty(varargin{1})
        adj = varargin{1};
    else
        adj = false;
    end
    
    if numel(varargin) >= 2 && ~isempty(varargin{2})
        tmpsz = varargin{2};
    else
        tmpsz = size(X(:,:,1)) + size(A(:,:,1)) - 1;
    end
    
    if numel(varargin) >= 3 && ~isempty(varargin{3})
        outsz = varargin{3};
    else
        outsz = size(X(:,:,1));
    end
    
    % Initialize output array
    Y = zeros([outsz, size(A,3)]);
    
    % Perform convolution for each slice
    for k = 1:size(A,3)
        if adj
            Y(:,:,k) = convfft2(A(:,:,k), X(:,:,k), true, tmpsz, outsz);
        else
            Y(:,:,k) = convfft2(A(:,:,k), X(:,:,k), false, tmpsz, outsz);
        end
    end
end 