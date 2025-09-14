function plot_2d_projections(metrics)
    % Define a linear colormap from #ee5253 to #10ac84
    num_datasets = metrics.num_datasets;
    cmap_start = [0.9333, 0.3216, 0.3255];  % RGB for #ee5253
    cmap_end = [0.0627, 0.6745, 0.5176];    % RGB for #10ac84
    cmap = [linspace(cmap_start(1), cmap_end(1), num_datasets)', ...
            linspace(cmap_start(2), cmap_end(2), num_datasets)', ...
            linspace(cmap_start(3), cmap_end(3), num_datasets)'];

    % Define marker shapes for different lambda or mini-loop values
    marker_shapes = {'o', 's', 'd', '^', 'v', '>', '<', 'p', 'h', '*', '+', 'x'};

    % Define control parameters for each dataset
    activation_levels = [3, 3, 1, 1, 1, 2, 2, 2, 2, 1, 1, 1];
    kernel_sizes = [3, 1, 3, 1, 1, 2, 2, 1, 1, 2, 2, 1];
    snr_levels = [1, 3, 3, 1, 3, 2, 3, 2, 3, 2, 3, 2];

    % Compute combined score if not already present
    if ~isfield(metrics, 'combined_score')
        metrics.combined_score = computeCombined_activationScore(metrics.demixing_score, ...
                                                metrics.activation_accuracy_final);
    end

    % Define a grayscale colormap for each control parameter
    cmap_activation = flipud(gray(max(activation_levels)));
    cmap_kernel = flipud(gray(max(kernel_sizes)));
    cmap_snr = flipud(gray(max(snr_levels)));

    % 1. Activation Levels
    fig_activation = figure('Position', [100 100 800 600], 'Name', 'Activation Levels');
    t_activation = tiledlayout(fig_activation, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    
    % Plot for activation levels
    ax1 = nexttile(t_activation);
    hold(ax1, 'on');
    for i = 1:num_datasets
        for j = 1:length(metrics.lambda1_values)
            plot(ax1, metrics.mini_loop_values, metrics.kernel_quality_final(i, :, j), 'o-', ...
                'Color', cmap_activation(activation_levels(i),:), 'Marker', marker_shapes{j});
        end
    end
    title(ax1, 'Kernel Similarity vs Mini-loop (Activation)');
    xlabel(ax1, 'Mini-loop');
    ylabel(ax1, 'Kernel Similarity');
    hold(ax1, 'off');

    ax2 = nexttile(t_activation);
    hold(ax2, 'on');
    for i = 1:num_datasets
        for j = 1:length(metrics.mini_loop_values)
            plot(ax2, metrics.lambda1_values, squeeze(metrics.kernel_quality_final(i, j, :)), 'o-', ...
                'Color', cmap_activation(activation_levels(i),:), 'Marker', marker_shapes{j});
        end
    end
    title(ax2, 'Kernel Similarity vs Lambda (Activation)');
    xlabel(ax2, 'Lambda');
    ylabel(ax2, 'Kernel Similarity');
    hold(ax2, 'off');

    ax3 = nexttile(t_activation);
    hold(ax3, 'on');
    for i = 1:num_datasets
        for j = 1:length(metrics.lambda1_values)
            plot(ax3, metrics.mini_loop_values, metrics.combined_score(i, :, j), 'o-', ...
                'Color', cmap_activation(activation_levels(i),:), 'Marker', marker_shapes{j});
        end
    end
    title(ax3, 'Combined Score vs Mini-loop (Activation)');
    xlabel(ax3, 'Mini-loop');
    ylabel(ax3, 'Combined Score');
    hold(ax3, 'off');

    ax4 = nexttile(t_activation);
    hold(ax4, 'on');
    for i = 1:num_datasets
        for j = 1:length(metrics.mini_loop_values)
            plot(ax4, metrics.lambda1_values, squeeze(metrics.combined_score(i, j, :)), 'o-', ...
                'Color', cmap_activation(activation_levels(i),:), 'Marker', marker_shapes{j});
        end
    end
    title(ax4, 'Combined Score vs Lambda (Activation)');
    xlabel(ax4, 'Lambda');
    ylabel(ax4, 'Combined Score');
    hold(ax4, 'off');

    % 2. Kernel Sizes
    fig_kernel = figure('Position', [100 100 800 600], 'Name', 'Kernel Sizes');
    t_kernel = tiledlayout(fig_kernel, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    
    % Plot for kernel sizes
    ax1 = nexttile(t_kernel);
    hold(ax1, 'on');
    for i = 1:num_datasets
        for j = 1:length(metrics.lambda1_values)
            plot(ax1, metrics.mini_loop_values, metrics.kernel_quality_final(i, :, j), 'o-', ...
                'Color', cmap_kernel(kernel_sizes(i),:), 'Marker', marker_shapes{j});
        end
    end
    title(ax1, 'Kernel Similarity vs Mini-loop (Kernel)');
    xlabel(ax1, 'Mini-loop');
    ylabel(ax1, 'Kernel Similarity');
    hold(ax1, 'off');

    ax2 = nexttile(t_kernel);
    hold(ax2, 'on');
    for i = 1:num_datasets
        for j = 1:length(metrics.mini_loop_values)
            plot(ax2, metrics.lambda1_values, squeeze(metrics.kernel_quality_final(i, j, :)), 'o-', ...
                'Color', cmap_kernel(kernel_sizes(i),:), 'Marker', marker_shapes{j});
        end
    end
    title(ax2, 'Kernel Similarity vs Lambda (Kernel)');
    xlabel(ax2, 'Lambda');
    ylabel(ax2, 'Kernel Similarity');
    hold(ax2, 'off');

    ax3 = nexttile(t_kernel);
    hold(ax3, 'on');
    for i = 1:num_datasets
        for j = 1:length(metrics.lambda1_values)
            plot(ax3, metrics.mini_loop_values, metrics.combined_score(i, :, j), 'o-', ...
                'Color', cmap_kernel(kernel_sizes(i),:), 'Marker', marker_shapes{j});
        end
    end
    title(ax3, 'Combined Score vs Mini-loop (Kernel)');
    xlabel(ax3, 'Mini-loop');
    ylabel(ax3, 'Combined Score');
    hold(ax3, 'off');

    ax4 = nexttile(t_kernel);
    hold(ax4, 'on');
    for i = 1:num_datasets
        for j = 1:length(metrics.mini_loop_values)
            plot(ax4, metrics.lambda1_values, squeeze(metrics.combined_score(i, j, :)), 'o-', ...
                'Color', cmap_kernel(kernel_sizes(i),:), 'Marker', marker_shapes{j});
        end
    end
    title(ax4, 'Combined Score vs Lambda (Kernel)');
    xlabel(ax4, 'Lambda');
    ylabel(ax4, 'Combined Score');
    hold(ax4, 'off');

    % 3. SNR Levels
    fig_snr = figure('Position', [100 100 800 600], 'Name', 'SNR Levels');
    t_snr = tiledlayout(fig_snr, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    
    % Plot for SNR levels
    ax1 = nexttile(t_snr);
    hold(ax1, 'on');
    for i = 1:num_datasets
        for j = 1:length(metrics.lambda1_values)
            plot(ax1, metrics.mini_loop_values, metrics.kernel_quality_final(i, :, j), 'o-', ...
                'Color', cmap_snr(snr_levels(i),:), 'Marker', marker_shapes{j});
        end
    end
    title(ax1, 'Kernel Similarity vs Mini-loop (SNR)');
    xlabel(ax1, 'Mini-loop');
    ylabel(ax1, 'Kernel Similarity');
    hold(ax1, 'off');

    ax2 = nexttile(t_snr);
    hold(ax2, 'on');
    for i = 1:num_datasets
        for j = 1:length(metrics.mini_loop_values)
            plot(ax2, metrics.lambda1_values, squeeze(metrics.kernel_quality_final(i, j, :)), 'o-', ...
                'Color', cmap_snr(snr_levels(i),:), 'Marker', marker_shapes{j});
        end
    end
    title(ax2, 'Kernel Similarity vs Lambda (SNR)');
    xlabel(ax2, 'Lambda');
    ylabel(ax2, 'Kernel Similarity');
    hold(ax2, 'off');

    ax3 = nexttile(t_snr);
    hold(ax3, 'on');
    for i = 1:num_datasets
        for j = 1:length(metrics.lambda1_values)
            plot(ax3, metrics.mini_loop_values, metrics.combined_score(i, :, j), 'o-', ...
                'Color', cmap_snr(snr_levels(i),:), 'Marker', marker_shapes{j});
        end
    end
    title(ax3, 'Combined Score vs Mini-loop (SNR)');
    xlabel(ax3, 'Mini-loop');
    ylabel(ax3, 'Combined Score');
    hold(ax3, 'off');

    ax4 = nexttile(t_snr);
    hold(ax4, 'on');
    for i = 1:num_datasets
        for j = 1:length(metrics.mini_loop_values)
            plot(ax4, metrics.lambda1_values, squeeze(metrics.combined_score(i, j, :)), 'o-', ...
                'Color', cmap_snr(snr_levels(i),:), 'Marker', marker_shapes{j});
        end
    end
    title(ax4, 'Combined Score vs Lambda (SNR)');
    xlabel(ax4, 'Lambda');
    ylabel(ax4, 'Combined Score');
    hold(ax4, 'off');

    % Add legends for marker shapes
    legend(ax2, arrayfun(@(x) sprintf('Mini-loop %d', x), metrics.mini_loop_values, 'UniformOutput', false));
    legend(ax4, arrayfun(@(x) sprintf('Mini-loop %d', x), metrics.mini_loop_values, 'UniformOutput', false));
end
