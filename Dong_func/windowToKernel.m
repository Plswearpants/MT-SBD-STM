function kernel_out = windowToKernel(kernel, window_type, varargin)
    % This function applies a window to a 2D square kernel
    % Inputs:
    %   kernel: 2D square matrix
    %   window_type: string, specifies the type of window
    %               'hann' - Hann window
    %               'hamming' - Hamming window
    %               'blackman' - Blackman window
    %               'gaussian' - Gaussian window (optional param: alpha)
    %               'kaiser' - Kaiser window (optional param: beta)
    %   varargin: optional parameters for specific windows
    %            gaussian: alpha (standard deviation, default = 2.5)
    %            kaiser: beta (shape parameter, default = 5)
    % Output:
    %   kernel_out: windowed kernel
    
    % Check if kernel is square
    [h, w] = size(kernel);
    if h ~= w
        error('Input kernel must be square');
    end
    
    % Create 2D window based on type
    switch lower(window_type)
        case 'hann'
            win_1d = hann(h);
            
        case 'hamming'
            win_1d = hamming(h);
            
        case 'blackman'
            win_1d = blackman(h);
            
        case 'gaussian'
            if ~isempty(varargin)
                alpha = varargin{1};  % custom alpha (standard deviation)
            else
                alpha = 2.5;  % default value
            end
            win_1d = gausswin(h, alpha);
            
        case 'kaiser'
            if ~isempty(varargin)
                beta = varargin{1};  % custom beta
            else
                beta = 5;  % default value
            end
            win_1d = kaiser(h, beta);
            
        otherwise
            error('Unsupported window type');
    end
    
    % Create 2D window by outer product
    window_2d = win_1d * win_1d';
    
    % Normalize window to preserve kernel energy
    window_2d = window_2d / sum(window_2d(:));
    
    % Apply window to kernel
    kernel_out = kernel .* window_2d;
    
    % Normalize output kernel to preserve energy
    kernel_out = kernel_out / sum(kernel_out(:)) * sum(kernel(:));
end
