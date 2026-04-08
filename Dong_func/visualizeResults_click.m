function Y_reconstructed = visualizeResults_click(Y, Aout, X0, Xout, bout, plot_title)
% Visualize clicked dataset results with a compact two-figure layout.
%
% Figure 1:
%   - Observation Y
%   - Total reconstructed observation sum_k conv(A_k, X_k) + sum_k b_k
%
% Figure 2:
%   - n x 3 grid for n kernels, per row:
%       [1] Per-kernel reconstruction: conv(A_k, X_k) + b_k
%       [2] Activation map X_k (Gaussian broadened)
%       [3] Kernel A_k
%
% Inputs:
%   Y          : observation image [H x W]
%   Aout       : reconstructed kernels (1 x n cell)
%   X0         : ground-truth activations [H x W x n] (used for alignment/filtering)
%   Xout       : reconstructed activations [H x W x n]
%   bout       : bias terms [n x 1] (or scalar/empty)
%   plot_title : optional title suffix

    if nargin < 5
        plot_title = '';
    end
    if isempty(Y) || isempty(Aout) || isempty(Xout)
        error('Y, Aout, and Xout are required.');
    end

    num_kernels = numel(Aout);
    if size(Xout, 3) ~= num_kernels
        error('Size mismatch: Xout channels must match number of kernels in Aout.');
    end

    if isempty(bout)
        bout = zeros(num_kernels, 1);
    else
        bout = bout(:);
        if numel(bout) == 1
            bout = repmat(bout, num_kernels, 1);
        elseif numel(bout) < num_kernels
            bout = [bout; zeros(num_kernels - numel(bout), 1)];
        elseif numel(bout) > num_kernels
            bout = bout(1:num_kernels);
        end
    end

    % Build kernel size table.
    kernel_sizes = zeros(num_kernels, 2);
    for k = 1:num_kernels
        kernel_sizes(k,:) = size(Aout{k});
    end

    % Match visualizeResults.m behavior:
    % align Xout to X0 first, then use filtered maps for activation display.
    use_alignment_pipeline = ~isempty(X0) && isequal(size(X0), size(Xout));
    X_for_recon = Xout;
    X_for_display = [];
    if use_alignment_pipeline
        try
            [~, aligned_maps] = evaluateActivationReconstruction(X0, Xout, kernel_sizes, false);
            if isfield(aligned_maps, 'Xout_aligned')
                X_for_recon = aligned_maps.Xout_aligned;
            end
            if isfield(aligned_maps, 'filtered') && ~isempty(aligned_maps.filtered)
                X_for_display = zeros(size(Xout));
                for k = 1:num_kernels
                    X_for_display(:,:,k) = aligned_maps.filtered(k).Xout;
                end
            end
        catch
            % Fallback below if alignment/similarity path fails for any reason.
            X_for_recon = Xout;
            X_for_display = [];
        end
    end

    % Fallback broadening if aligned filtered maps are unavailable.
    if isempty(X_for_display)
        X_for_display = gaussian_broaden_channels(X_for_recon, kernel_sizes);
    end

    % Build per-kernel and total reconstructions.
    Y_per_kernel = zeros([size(Y), num_kernels]);
    Y_reconstructed = zeros(size(Y));
    for k = 1:num_kernels
        Yk = convfft2(Aout{k}, X_for_recon(:,:,k)) + bout(k);
        Y_per_kernel(:,:,k) = Yk;
        Y_reconstructed = Y_reconstructed + Yk;
    end

    % Shared display limits for Y-like images.
    y_stack = cat(3, Y, Y_reconstructed, Y_per_kernel);
    y_min = min(y_stack(:));
    y_max = max(y_stack(:));
    if y_min == y_max
        y_max = y_min + 1;
    end

    % Figure 1: observation vs total reconstruction.
    figure('Name', 'Observation vs Reconstructed', 'Position', [80 80 980 430]);
    subplot(1,2,1);
    imagesc(Y);
    axis image off;
    caxis([y_min y_max]);
    colormap(gca, gray);
    colorbar;
    title('Observation');

    subplot(1,2,2);
    imagesc(Y_reconstructed);
    axis image off;
    caxis([y_min y_max]);
    colormap(gca, gray);
    colorbar;
    title('Reconstructed Observation');

    if ~isempty(plot_title)
        sgtitle(sprintf('Observation Pair | %s', plot_title), 'FontWeight', 'bold');
    else
        sgtitle('Observation Pair', 'FontWeight', 'bold');
    end

    % Figure 2: n x 3 grid (Y_k, X_k broadened, A_k).
    figure('Name', 'Per-Kernel Decomposition', ...
           'Position', [120 80 1200 max(320, 260 * num_kernels)]);

    for k = 1:num_kernels
        % Column 1: per-kernel reconstructed observation
        subplot(num_kernels, 3, (k-1)*3 + 1);
        imagesc(Y_per_kernel(:,:,k));
        axis image off;
        caxis([y_min y_max]);
        colormap(gca, gray);
        colorbar;
        title(sprintf('Y_{rec}^{(%d)} = A_{%d} * X_{%d} + b_{%d}', k, k, k, k), ...
            'Interpreter', 'tex');

        % Column 2: broadened activation map
        subplot(num_kernels, 3, (k-1)*3 + 2);
        imagesc(X_for_display(:,:,k));
        axis image off;
        colormap(gca, hot);
        colorbar;
        title(sprintf('X_{%d} (filtered for display)', k), 'Interpreter', 'tex');

        % Column 3: kernel
        subplot(num_kernels, 3, (k-1)*3 + 3);
        imagesc(Aout{k});
        axis image off;
        colormap(gca, gray);
        colorbar;
        title(sprintf('A_{%d}', k), 'Interpreter', 'tex');
    end

    if ~isempty(plot_title)
        sgtitle(sprintf('Per-Kernel Components | %s', plot_title), 'FontWeight', 'bold');
    else
        sgtitle('Per-Kernel Components', 'FontWeight', 'bold');
    end
end

function Xout_gaussian = gaussian_broaden_channels(X_out, kernel_sizes)
% Apply per-channel Gaussian broadening using kernel-size-dependent sigma.
    [h, w, num_kernels] = size(X_out);
    Xout_gaussian = zeros(h, w, num_kernels);
    for k = 1:num_kernels
        sigma = min(kernel_sizes(k,:)) / 10;
        sigma = max(sigma, 1e-6);
        window_size = max(1, ceil(3 * sigma));
        [x, y] = meshgrid(-window_size:window_size);
        gaussian_kernel = exp(-(x.^2 + y.^2) / (2 * sigma^2));
        gaussian_kernel = gaussian_kernel / sum(gaussian_kernel(:));
        Xout_gaussian(:,:,k) = conv2(X_out(:,:,k), gaussian_kernel, 'same');
    end
end
