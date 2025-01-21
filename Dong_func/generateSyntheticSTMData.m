function [Y, A0_noiseless, A0, X0, SNR, params] = generateSyntheticSTMData(varargin)
    % Parse required input parameters
    p = inputParser;
    p.addRequired('num_kernels');
    p.addRequired('image_size');
    p.addRequired('kernel_sizes');
    p.addRequired('SNR');
    p.addRequired('theta_cap');
    p.addRequired('interactive_selection');
    p.addRequired('LDoS_path');
    
    % Parse inputs
    p.parse(varargin{:});
    
    % Extract parameters
    num_kernels = p.Results.num_kernels;
    image_size = p.Results.image_size;
    kernel_sizes = p.Results.kernel_sizes;
    SNR = p.Results.SNR;
    theta_cap = p.Results.theta_cap;
    interactive_selection = p.Results.interactive_selection;
    LDoS_path = p.Results.LDoS_path;
    
    % Validate inputs
    validateattributes(num_kernels, {'numeric'}, {'positive', 'integer', 'scalar'});
    validateattributes(image_size, {'numeric'}, {'positive', 'integer', 'size', [1, 2]});
    validateattributes(kernel_sizes, {'numeric'}, {'positive', 'integer', 'size', [num_kernels, 2]});
    validateattributes(SNR, {'numeric'}, {'positive', 'scalar'});
    validateattributes(theta_cap, {'numeric'}, {'positive', 'scalar'});
    validateattributes(interactive_selection, {'logical'}, {'scalar'});
    
    % Load LDoS data
    if ~exist(LDoS_path, 'file')
        error('LDoS data file not found: %s', LDoS_path);
    end
    load(LDoS_path);
    
    % Get kernels either interactively or randomly
    if interactive_selection
        [A0, A0_noiseless, kernel_params] = selectKernelsInteractive(LDoS_sim);
        % Update parameters based on user selection
        num_kernels = kernel_params.num_kernels;
        kernel_sizes = kernel_params.kernel_size;
        image_size = kernel_params.image_size;
        sliceidx = kernel_params.sliceidx;
        eta_kernel = kernel_params.eta_kernel;
    else
        % Random selection process
        sliceidx = randperm(size(LDoS_sim,3), num_kernels);
        A0 = cell(1,num_kernels);
        A0_noiseless = cell(1,num_kernels);
        
        % Calculate average variance and prepare noiseless kernels
        avg_var = 0;
        for n = 1:num_kernels
            A0_noiseless{n} = imresize(LDoS_sim(:,:,sliceidx(n)), kernel_sizes(n,:));
            A0_noiseless{n} = proj2oblique(A0_noiseless{n});
            avg_var = avg_var + var(A0_noiseless{n},0,"all");
        end
        avg_var = avg_var / num_kernels;
        
        % Apply universal noise to kernels
        eta_kernel = avg_var/SNR;
        for n = 1:num_kernels
            A0{n} = A0_noiseless{n} + sqrt(eta_kernel)*randn(size(A0_noiseless{n}));
            A0{n} = proj2oblique(A0{n});
        end
    end
    
    % Generate activation maps
    theta = theta_cap/2 + theta_cap/2 * rand(1, num_kernels);
    X0 = zeros([image_size num_kernels]);
    for k = 1:num_kernels
        X0_good = false;
        while ~X0_good
            X0(:,:,k) = double(rand(image_size) <= theta(k));
            X0_good = sum(X0(:,:,k) ~= 0) > 0;
        end
    end
    
    % Generate clean observation and add noise
    Y_clean = zeros(size(X0, 1:2));
    for k = 1:size(X0, 3)
        Y_clean = Y_clean + convfft2(A0_noiseless{k}, X0(:,:,k));
    end
    
    eta = var(Y_clean, 0, "all") / SNR;
    Y = Y_clean + sqrt(eta) * randn(size(Y_clean));
    
    % Package parameters for output
    params = struct();
    params.num_kernels = num_kernels;
    params.image_size = image_size;
    params.kernel_sizes = kernel_sizes;
    params.SNR = SNR;
    params.theta_cap = theta_cap;
    params.theta = theta;
    params.eta = eta;
    params.eta_kernel = eta_kernel;
    params.selected_slices = sliceidx;
    if interactive_selection
        params.kernel_params = kernel_params;
    end
end