function plot_defects_snr_kernel_similarity(dataset_metrics)
% Plots a 2D scatter plot: x = number of defects, y = kernel similarity, color = SNR
% dataset_metrics: structure with fields X0, SNR_values, kernel_quality_final

% Prepare arrays to collect data
num_points = numel(dataset_metrics.kernel_quality_final);
defects = nan(num_points, 1);
SNRs = nan(num_points, 1);
kernel_sim = nan(num_points, 1);

% Loop through all parameter combinations
count = 0;
for i = 1:numel(dataset_metrics.kernel_quality_final)
    kq = dataset_metrics.kernel_quality_final(i);
    if isnan(kq)
        continue;
    end
    % Find the corresponding X0 and SNR
    [idx1, idx2, idx3] = ind2sub(size(dataset_metrics.kernel_quality_final), i);
    X0 = dataset_metrics.X0{idx1, idx2, idx3};
    if isempty(X0)
        continue;
    end
    % Compute average number of defects (sum of each activation map, then mean over kernels)
    num_kernels = size(X0, 3);
    norms = zeros(1, num_kernels);
    for k = 1:num_kernels
        norms(k) = sum(X0(:,:,k), 'all');
    end
    avg_defects = mean(norms);
    % Get SNR value
    SNR = dataset_metrics.SNR_values(idx1);
    % Store
    count = count + 1;
    defects(count) = avg_defects;
    SNRs(count) = SNR;
    kernel_sim(count) = kq;
end

defects = defects(1:count);
SNRs = SNRs(1:count);
kernel_sim = kernel_sim(1:count);

% Create 2D scatter plot with SNR color coding
figure;
scatter(defects, kernel_sim, 60, SNRs, 'filled');
colormap(flipud(gray)); % dark = low SNR, light = high SNR
colorbar;
c = colorbar;
c.Label.String = 'SNR';
xlabel('Average Number of Defects');
xscale('log');
ylabel('Kernel Similarity');
title('Defects vs Kernel Similarity (SNR color-coded)');
grid on;
ylim([0 1]);
end 