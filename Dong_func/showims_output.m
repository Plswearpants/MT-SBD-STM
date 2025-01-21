function output_image = showims_output(Y, A0, X0, A, X, k, kplus, idx)
    %SHOWIMS_OUTPUT    Show images after each iteration and output as an image.
    A = reshape(A, [k size(Y,3)]);

    Y_hat = convfft2(A(:,:,idx),X);
    
    % Create a new figure for capturing
    fig = figure;
    
    % Set the figure size for higher resolution
    fig.Position = [100, 100, 1200, 800];  % [left, bottom, width, height]
    
    % Plot the images in subplots
    subplot(321); imagesc((Y(:,:,idx)));     ylabel('Y'); title('Original');
    subplot(322); imagesc((Y_hat));          title('Recovered');
    
    subplot(324); imagesc(abs(A(:,:,idx)));
    subplot(323); imagesc(abs(A0(:,:,idx)));
    
    if ~isempty(kplus)      % i.e. we're in Phase II
        ylabel('abs(A)  (lifted)'); 
        X_hat = circshift(X,kplus);
    else
        ylabel('abs(A)'); 
        X_hat = X;
    end
    
    subplot(325); imagesc(abs(X0)); ylabel('abs(X)');
    subplot(326); imagesc(abs(X_hat));
    
    % Adjust colormap and axis properties for better visualization
    colormap(gray);
    axis image;
    
    % Ensure the figure is rendered at the specified resolution
    drawnow;

    % Capture the figure content as a frame
    frame = getframe(fig);
    
    % Convert the frame to an image
    output_image = frame2im(frame);
    
    % Close the figure after capturing
    close(fig);
end
