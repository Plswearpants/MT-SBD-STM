function showims(Y, A0, X0, A, X, k, kplus, idx)
%SHOWIMS    Show images after each iteration, with square subplots and minimal margins.
%   Optimized for speed using persistent handles, tiledlayout, and direct CData updates.
%   Requires MATLAB R2019b+ for tiledlayout.

    persistent fig_handle img_handles
    
    % Compute data
    A = reshape(A, [k size(Y,3)]);
    Y_hat = convfft2(A(:,:,idx), X);
    if ~isempty(kplus)
        X_hat = circshift(X, kplus);
    else
        X_hat = X;
    end
    
    % First call: create figure and layout once
    if isempty(fig_handle) || ~isvalid(fig_handle)
        fig_handle = gcf;
        clf(fig_handle);
        set(fig_handle, 'Color', 'w', 'Renderer', 'painters');
        
        % Create tiled layout (faster than subplot/axes)
        tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
        
        % Pre-create all image objects with initial data
        % Row 1: Original observation | Initial guess kernel | GT activation
        nexttile;
        img_handles(1) = imagesc(Y(:,:,idx));
        axis image off; title('Y (Original)');
        
        nexttile;
        img_handles(2) = imagesc(abs(A0(:,:,idx)));
        axis image off; title('A₁ (Initial Guess)');
        
        nexttile;
        img_handles(3) = imagesc(abs(X0));
        axis image off; title('X₀ (GT Activation)');
        
        % Row 2: Recovered observation | Current estimated kernel | Current estimated activation
        nexttile;
        img_handles(4) = imagesc(Y_hat);
        axis image off; title('Ŷ_channel (Recovered)');
        
        nexttile;
        img_handles(5) = imagesc(abs(A(:,:,idx)));
        axis image off;
        if ~isempty(kplus)
            title('A (Est. Kernel, lifted)');
        else
            title('A (Est. Kernel)');
        end
        
        nexttile;
        img_handles(6) = imagesc(abs(X_hat));
        axis image off; title('X̂ (Est. Activation)');
        
        colormap(fig_handle, "gray");
    else
        % Subsequent calls: only update image data (fast!)
        set(img_handles(1), 'CData', Y(:,:,idx));
        set(img_handles(2), 'CData', abs(A0(:,:,idx)));  % Initial guess (stays constant)
        set(img_handles(3), 'CData', abs(X0));           % GT activation (stays constant)
        set(img_handles(4), 'CData', Y_hat);              % Recovered observation (updates)
        set(img_handles(5), 'CData', abs(A(:,:,idx)));    % Current estimated kernel (updates)
        set(img_handles(6), 'CData', abs(X_hat));         % Current estimated activation (updates)
    end
    
    drawnow limitrate;  % Fast update with rate limiting
end
