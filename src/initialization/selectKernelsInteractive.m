function [A0, A0_noiseless, kernel_params] = selectKernelsInteractive(LDoS_sim)
% SELECTKERNELSINTERACTIVE Interactive kernel selection and parameter definition
%   
%   Input:
%       LDoS_sim: 3D array of LDoS simulation data
%   
%   Outputs:
%       A0: Array of selected kernels with noise
%       A0_noiseless: Array of noiseless kernels
%       kernel_params: Struct containing kernel-related parameters

    % Display the 3D LDoS data for selection
    fprintf('Displaying 3D LDoS simulation data...\n');
    d3gridDisplay(LDoS_sim, 'dynamic');
    title('LDoS Simulation Data - Use for Kernel Selection');

    % Get user input for kernel selection
    num_kernels = input('\nEnter number of kernels (default=3): ');
    if isempty(num_kernels)
        num_kernels = 3;
    end

    % Select slices
    fprintf('\nPlease select %d slice indices for kernels (1-%d):\n', num_kernels, size(LDoS_sim,3));
    sliceidx = zeros(1, num_kernels);
    for k = 1:num_kernels
        sliceidx(k) = input(sprintf('Enter slice index for kernel %d: ', k));
        
        % Validate input
        while sliceidx(k) < 1 || sliceidx(k) > size(LDoS_sim,3)
            fprintf('Invalid slice index. Please enter a number between 1 and %d\n', size(LDoS_sim,3));
            sliceidx(k) = input(sprintf('Enter slice index for kernel %d: ', k));
        end
    end

    % Display selected slices for confirmation
    figure('Name', 'Selected Kernel Slices');
    for k = 1:num_kernels
        subplot(1, num_kernels, k);
        imagesc(LDoS_sim(:,:,sliceidx(k)));
        title(sprintf('Kernel %d (Slice %d)', k, sliceidx(k)));
        colorbar;
        axis square;
    end

    % Confirm selection
    confirmation = input('\nAre you satisfied with these kernel selections? (y/n): ', 's');
    while ~strcmpi(confirmation, 'y')
        % If not satisfied, allow reselection
        fprintf('\nPlease reselect slices:\n');
        for k = 1:num_kernels
            sliceidx(k) = input(sprintf('Enter slice index for kernel %d: ', k));
            while sliceidx(k) < 1 || sliceidx(k) > size(LDoS_sim,3)
                fprintf('Invalid slice index. Please enter a number between 1 and %d\n', size(LDoS_sim,3));
                sliceidx(k) = input(sprintf('Enter slice index for kernel %d: ', k));
            end
        end
        
        % Update display
        for k = 1:num_kernels
            subplot(1, num_kernels, k);
            imagesc(LDoS_sim(:,:,sliceidx(k)));
            title(sprintf('Kernel %d (Slice %d)', k, sliceidx(k)));
            colorbar;
            axis square;
        end
        
        confirmation = input('\nAre you satisfied with these kernel selections? (y/n): ', 's');
    end

    % Get additional parameters
    fprintf('\nDefining parameters...\n');

    % Image size
    image_size = input('Enter image size [width height] (default=[300,300]): ');
    if isempty(image_size)
        image_size = [300, 300];
    end

    % Kernel sizes
    kernel_size = zeros(num_kernels,2);
    for k = 1:num_kernels
        fprintf('\nKernel %d size:\n', k);
        size_input = input(sprintf('Enter size [width height] (default=[70,70]): '));
        if isempty(size_input)
            kernel_size(k,:) = [70, 70];
        else
            kernel_size(k,:) = size_input;
        end
    end

    % Initialize kernel arrays
    A0 = cell(1,num_kernels);
    A0_noiseless = cell(1,num_kernels);
    
    % Calculate average variance and prepare noiseless kernels
    avg_var = 0;
    for n = 1:num_kernels
        A0_noiseless{n} = imresize(LDoS_sim(:,:,sliceidx(n)), kernel_size(n,:));
        A0_noiseless{n} = proj2oblique(A0_noiseless{n});
        avg_var = avg_var + var(A0_noiseless{n},0,'all');
    end
    avg_var = avg_var / num_kernels;
    
    % Get SNR for noise addition
    SNR = input('\nEnter SNR for kernel noise (default=5): ');
    if isempty(SNR)
        SNR = 5;
    end
    
    % Apply universal noise to kernels
    eta_kernel = avg_var/SNR;
    for n = 1:num_kernels
        A0{n} = A0_noiseless{n} + sqrt(eta_kernel)*randn(size(A0_noiseless{n}));
        A0{n} = proj2oblique(A0{n});
    end

    % Package all kernel-related parameters
    kernel_params = struct('num_kernels', num_kernels, 'sliceidx', sliceidx, 'kernel_size', kernel_size, ...
        'image_size', image_size,'SNR', SNR, 'eta_kernel', eta_kernel); 
end