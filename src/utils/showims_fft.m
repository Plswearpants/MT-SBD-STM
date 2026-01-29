function showims_fft(Y, A0, X0, A, X, k, kplus, idx, saveImage)
    %SHOWIMS_FFT Show images after each iteration, with options for display size and saving.
    A = reshape(A, [k size(Y, 3)]);

    % Compute the recovered Y
    Y_hat = convfft2(A(:,:,idx), X);

    % Create a figure with a larger display size
    figure('Position', [0, 0, 1200, 1200]);

    % Ensure each subplot is square
    subplot(3, 3, 1); imagesc(Y(:,:,idx)); ylabel('Y'); title('original Y'); axis square;
    subplot(3, 3, 2); imagesc(Y_hat); title('Recovered Y'); axis square;
    subplot(3, 3, 3); imagesc(abs(fftshift(fft2(Y(:,:,idx) - mean(mean(Y(:,:,idx))))))); title('original Y FFT'); axis square;

    subplot(3, 3, 4); imagesc(abs(A0(:,:,idx))); title('Initial Kernel'); axis square;
    subplot(3, 3, 5); imagesc(abs(A(:,:,idx))); title('recovered Kernel'); axis square;
    subplot(3, 3, 6); imagesc(abs(fftshift(fft2(A(:,:,idx) - mean(mean(A(:,:,idx))))))); title('recovered Kernel FFT'); axis square;

    if ~isempty(kplus)  % i.e., we're in Phase II
        ylabel('abs(A) (lifted)'); 
        X_hat = circshift(X, kplus);
    else
        ylabel('abs(A)'); 
        X_hat = X;
    end

    subplot(3, 3, 7); imagesc(abs(X0)); ylabel('abs(X)'); title('recovered activation map'); axis square;
    subplot(3, 3, 8); imagesc(abs(X_hat)); title('recovered activation map'); axis square;
    subplot(3, 3, 9); imagesc(abs(fftshift(fft2(A0(:,:,idx) - mean(mean(A0(:,:,idx))))))); title('initial Kernel FFT'); axis square;

    % Save the figure if saveImage is true
    if saveImage
        filename = sprintf('showims_fft_output_%d.png', idx);
        saveas(gcf, filename);
    end

    drawnow;
end
