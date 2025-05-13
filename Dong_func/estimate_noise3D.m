function noise_variance = estimate_noise3D(image, method)
    % Accepts 2D or 3D image data. For 3D, applies the same ROI to all slices.
    if ndims(image) == 3
        num_slices = size(image, 3);
    else
        num_slices = 1;
    end

    % Use the first slice for ROI selection
    if num_slices > 1
        img_gray = image(:,:,1);
    else
        img_gray = image;
    end

    % Check the selected method and proceed accordingly
    switch lower(method)
        case 'std'
            % Display the first slice for ROI selection
            figure;
            imshow(img_gray, []);
            title('Draw a rectangle to select the region of interest (ROI) of noise. Double-click to confirm.');

            % Let the user draw a rectangle to select the ROI
            h = imrect;
            position = wait(h);  % Wait until the ROI is double-clicked

            % Preallocate output
            noise_variance = zeros(1, num_slices);

            % Apply the ROI to all slices
            for i = 1:num_slices
                if num_slices > 1
                    current_slice = image(:,:,i);
                else
                    current_slice = image;
                end
                roi = imcrop(current_slice, position);
                noise_variance(i) = var(double(roi(:)));
            end

            % Display the result
            %disp(['Estimated noise variance (standard deviation method) for each slice: ', mat2str(noise_variance)]);

            % Close the figure
            close;

        case 'wavelet'
            % Preallocate output
            noise_variance = zeros(1, num_slices);

            for i = 1:num_slices
                if num_slices > 1
                    img_gray = image(:,:,i);
                else
                    img_gray = image;
                end
                img_double = double(img_gray);

                % Perform wavelet decomposition using the 'db1' wavelet (Haar wavelet)
                [~, cH, cV, cD] = dwt2(img_double, 'db1');

                % Calculate the variance of the high-frequency components
                noise_var_H = var(cH(:));
                noise_var_V = var(cV(:));
                noise_var_D = var(cD(:));

                % Combine the variances as a measure of the noise level
                noise_variance(i) = noise_var_H + noise_var_V + noise_var_D;
            end

            % Display the result
            disp(['Estimated noise variance (wavelet method) for each slice: ', mat2str(noise_variance)]);

        otherwise
            error('Invalid method. Choose either \"std\" or \"wavelet\".');
    end
end
