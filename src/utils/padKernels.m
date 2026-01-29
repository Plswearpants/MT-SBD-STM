function A0_padded = padKernels(A0_noiseless, SNR, target_sizes)
% PADKERNELS Pads or crops kernels to target sizes and adds noise based on SNR
%   A0_noiseless: num_slices x num_kernel struct containing kernels
%   SNR: Signal-to-noise ratio for noise addition
%   target_sizes: num_kernel x 2 matrix of target [x,y] dimensions
%
%   Returns A0_padded as num_kernel x 1 struct, where each entry is a
%   target_size x num_slices matrix

% Get dimensions
num_slices = size(A0_noiseless, 1);
num_kernels = size(A0_noiseless, 2);

% Initialize output struct
A0_padded = cell(num_kernels,1);

% Process each kernel
for kernel_idx = 1:num_kernels
    % Get target size for this kernel
    target_size = target_sizes(kernel_idx, :);
    % Initialize 3D matrix for this kernel (target_size x num_slices)
    kernel_3D = zeros([target_size(1),target_size(2), num_slices]);
    % Process each slice
    for slice_idx = 1:num_slices
        % Get current kernel
        current_kernel = A0_noiseless{slice_idx, kernel_idx};
        current_size = size(current_kernel);
        
        if all(target_size >= current_size)
            % Case 1: Target size is larger - pad with zeros
            % Calculate padding amounts
            pad_x = floor((target_size(1) - current_size(1)) / 2);
            pad_y = floor((target_size(2) - current_size(2)) / 2);
            
            % Place original kernel in center
            kernel_3D(pad_x+1:pad_x+current_size(1), ...
                     pad_y+1:pad_y+current_size(2), ...
                     slice_idx) = current_kernel;
            
            % Calculate signal power (variance of the original kernel)
            signal_power = var(current_kernel(:));
            
        else
            % Case 2: Target size is smaller - crop from center
            % Calculate crop amounts
            crop_x = floor((current_size(1) - target_size(1)) / 2);
            crop_y = floor((current_size(2) - target_size(2)) / 2);
            
            % Extract center portion
            kernel_3D(:,:,slice_idx) = current_kernel(crop_x+1:crop_x+target_size(1), ...
                                                    crop_y+1:crop_y+target_size(2));
            
            % Calculate signal power (variance of the cropped kernel)
            current_kernel = kernel_3D(:,:,slice_idx);
            signal_power = var(current_kernel(:));
        end
        
        % Calculate noise power based on SNR
        noise_power = signal_power / SNR;
        
        % Generate noise with calculated power
        noise = sqrt(noise_power) * randn(target_size);
        
        % Add noise to padded/cropped kernel
        kernel_3D(:,:,slice_idx) = kernel_3D(:,:,slice_idx) + noise;
    end
    
    % Store in output struct
    A0_padded{kernel_idx} = kernel_3D;
end
end 