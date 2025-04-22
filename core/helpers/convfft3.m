function [ Y ] = convfft3( A, X, varargin )
%CONVFFT3   FFT implementation of 2D linear convolution for 3D arrays
%   Performs 2D convolution along the third dimension of 3D arrays
%
%   Y = convfft3(A, X) convolves each slice of A with corresponding slice of X
%   A and X must be 3D arrays with matching third dimension
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
    if ndims(A) ~= 3 || ndims(X) ~= 3
        error('Inputs must be 3D arrays');
    end
    
    if size(A,3) ~= size(X,3)
        error('Third dimension of A and X must match');
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