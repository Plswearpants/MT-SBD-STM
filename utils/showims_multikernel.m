function showims_multikernel(Y, A0, X0, A, X, k, kplus, idx)
    %SHOWIMS_MULTIKERNEL     Show images after each iteration for multi-kernel case.
    A = reshape(A, [k size(A,3)]);
    num_kernels = size(A, 3);
    
    % Calculate Y_hat using all kernels
    Y_hat = zeros(size(Y));
    for i = 1:num_kernels
        Y_hat = Y_hat + convfft2(A(:,:,i), X(:,:,i));
    end
    
    % Plot kernels (A0 and A)
    num_rows = num_kernels+1;
    figure;

    % Plot Y and Y_hat
    subplot(num_rows, 4, 1); imagesc(Y); ylabel('Y'); title('Original');
    subplot(num_rows, 4, 3); imagesc(Y_hat); title('Recovered');

    if ~isempty(kplus)      % i.e. we're in Phase II
        X_hat = circshift(X(:,:,i), kplus);
    else
        X_hat = X(:,:,i);
    end

    for i = 1:num_kernels
        subplot(num_rows, 4, 4*i+1); 
        imagesc(abs(A0(:,:,i))); 
        title(sprintf('Original Kernel %d', i));
        subplot(num_rows, 4, 4*i+2);
        imagesc(abs(X0(:,:,i)));
        title(sprintf('Original X %d', i));
        subplot(num_rows, 4, 4*i+3); 
        imagesc(abs(A(:,:,i))); 
        title(sprintf('Recovered Kernel %d', i));
        subplot(num_rows, 4, 4*i+4);
        imagesc(abs(X_hat));
        title(sprintf('Recovered X %d', i));
    end
    
    drawnow;
end