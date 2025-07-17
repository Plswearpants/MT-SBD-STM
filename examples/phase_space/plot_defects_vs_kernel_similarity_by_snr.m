function plot_defects_vs_kernel_similarity_by_snr(combined_metrics)
% For each unique SNR, plot a 2D scatter: x = avg number of defects, y = kernel similarity

SNRs = combined_metrics.SNR_values;
num_SNR = length(SNRs);

% Prepare data
num_points = numel(combined_metrics.kernel_quality_final);
defects = nan(num_points, 1);
kernel_sim = nan(num_points, 1);
SNR_idx = nan(num_points, 1);

count = 0;
for i = 1:numel(combined_metrics.kernel_quality_final)
    kq = combined_metrics.kernel_quality_final(i);
    if isnan(kq)
        continue;
    end
    [idx1, idx2, idx3] = ind2sub(size(combined_metrics.kernel_quality_final), i);
    X0 = combined_metrics.X0{idx1, idx2, idx3};
    if isempty(X0)
        continue;
    end
    num_kernels = size(X0, 3);
    norms = zeros(1, num_kernels);
    for k = 1:num_kernels
        norms(k) = sum(X0(:,:,k), 'all');
    end
    avg_defects = mean(norms);
    count = count + 1;
    defects(count) = avg_defects;
    kernel_sim(count) = kq;
    SNR_idx(count) = idx1;
end

defects = defects(1:count);
kernel_sim = kernel_sim(1:count);
SNR_idx = SNR_idx(1:count);

% Only keep points with kernel similarity > 0.84
valid = kernel_sim > 0.84;
defects = defects(valid);
kernel_sim = kernel_sim(valid);
SNR_idx = SNR_idx(valid);

% After preparing defects, kernel_sim, SNR_idx, etc.
x_min = min(defects);
x_max = max(defects);
cmap_data = load('metric_colormapv3.mat');


% Plot: horizontal subplots, sharing the same Y axis
figure;
t = tiledlayout(1, num_SNR, 'TileSpacing', 'compact', 'Padding', 'compact');
for s = 1:num_SNR
    colormap(cmap_data.CustomColormap);
    ax = nexttile;
    idx = SNR_idx == s;
    % Plot all points with color according to kernel similarity
    scatter(defects(idx), kernel_sim(idx), 60, kernel_sim(idx), 'filled');
    c = colorbar;
    clim([0 1]);
    c.Label.String = 'Kernel Similarity';
    xlabel('Average Number of Defects');
    xscale('log');
    ylabel('Kernel Similarity');
    title(sprintf('SNR = %.2f', SNRs(s)));
    grid on;
    ylim([0.84 1]);
    xlim([x_min, x_max]);
    if s > 1
        ax.YTickLabel = [];
        ylabel('');
    end
end
sgtitle('Defects vs Kernel Similarity by SNR');
end 