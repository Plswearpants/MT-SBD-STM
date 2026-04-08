function metrics2heat_multiple_runs(dataset_metrics_array, run_names, mode)
    % Create 3D heatspace visualization for multiple experiment runs.
    % Inputs:
    %   dataset_metrics_array : cell array of dataset_metrics structures
    %   run_names             : cell array of names for each run (optional)
    %   mode                  : 1 -> axis-3 is Nobs (default), 2 -> side_length_ratio
    
    if nargin < 2 || isempty(run_names)
        run_names = cell(length(dataset_metrics_array), 1);
        for i = 1:length(dataset_metrics_array)
            run_names{i} = sprintf('Run %d', i);
        end
    end
    if nargin < 3 || isempty(mode)
        mode = 1;
    end
    
    if length(dataset_metrics_array) ~= length(run_names)
        error('Number of dataset_metrics must match number of run_names');
    end
    
    % Select metric field names based on mode
    if mode == 2
        kernel_field = 'kernel_quality_final_by_side_length_ratio';
        combined_field = 'combined_activationScore_by_side_length_ratio';
        axis3_source = 'side_length_ratio_values';
        axis3_label = 'side length ratio';
        axis3_tick_fmt = '%.3f';
    else
        kernel_field = 'kernel_quality_final';
        combined_field = 'combined_activationScore';
        axis3_source = 'Nobs_values';
        axis3_label = 'N_{obs}';
        axis3_tick_fmt = '%.0f';
    end
    
    figure('Position', [100 100 1600 800]);
    
    subplot(1,2,1);
    plot_3D_heatspace_multiple(dataset_metrics_array, run_names, ...
        kernel_field, 'Kernel Similarity', axis3_source, axis3_label, axis3_tick_fmt);
    
    subplot(1,2,2);
    plot_3D_heatspace_multiple(dataset_metrics_array, run_names, ...
        combined_field, 'Combined Activation Score', axis3_source, axis3_label, axis3_tick_fmt);
end

function h = plot_3D_heatspace_multiple(dataset_metrics_array, run_names, metric_field, title_str, axis3_source, axis3_label, axis3_tick_fmt)
    
    all_snr_values = [];
    all_theta_values = [];
    all_axis3_values = [];
    
    for i = 1:length(dataset_metrics_array)
        metrics = dataset_metrics_array{i};
        all_snr_values = [all_snr_values; metrics.SNR_values(:)];
        all_theta_values = [all_theta_values; metrics.theta_cap_values(:)];
        if isfield(metrics, axis3_source)
            all_axis3_values = [all_axis3_values; metrics.(axis3_source)(:)];
        end
    end
    
    unique_snr = unique(all_snr_values);
    unique_theta = unique(all_theta_values);
    unique_axis3 = unique(all_axis3_values);
    
    [SNR, Theta, Axis3] = meshgrid(unique_snr, unique_theta, unique_axis3);
    
    SNR = permute(SNR, [2, 1, 3]);
    Theta = permute(Theta, [2, 1, 3]);
    Axis3 = permute(Axis3, [2, 1, 3]);
    
    num_snr_values = length(unique_snr);
    if num_snr_values > 1
        snr_spaced = linspace(min(unique_snr), max(unique_snr), num_snr_values);
        snr_mapping = containers.Map(unique_snr, snr_spaced);
        SNR_remapped = zeros(size(SNR));
        for i = 1:numel(SNR)
            SNR_remapped(i) = snr_mapping(SNR(i));
        end
        SNR = SNR_remapped;
    end
    
    combined_metric = nan(size(SNR));
    point_counts = zeros(size(SNR));
    
    for run_idx = 1:length(dataset_metrics_array)
        metrics = dataset_metrics_array{run_idx};
        
        if ~isfield(metrics, metric_field)
            warning('Field %s not found in run %d', metric_field, run_idx);
            continue;
        end
        run_metric = average_over_reps(metrics.(metric_field));
        
        run_axis3 = metrics.(axis3_source);
        
        for i = 1:length(metrics.SNR_values)
            for j = 1:length(metrics.theta_cap_values)
                for k = 1:length(run_axis3)
                    [~, snr_idx] = min(abs(unique_snr - metrics.SNR_values(i)));
                    [~, theta_idx] = min(abs(unique_theta - metrics.theta_cap_values(j)));
                    [~, axis3_idx] = min(abs(unique_axis3 - run_axis3(k)));
                    
                    val = run_metric(i, j, k);
                    if ~isnan(val)
                        if isnan(combined_metric(snr_idx, theta_idx, axis3_idx))
                            combined_metric(snr_idx, theta_idx, axis3_idx) = val;
                            point_counts(snr_idx, theta_idx, axis3_idx) = 1;
                        else
                            n = point_counts(snr_idx, theta_idx, axis3_idx);
                            combined_metric(snr_idx, theta_idx, axis3_idx) = ...
                                (combined_metric(snr_idx, theta_idx, axis3_idx) * n + val) / (n + 1);
                            point_counts(snr_idx, theta_idx, axis3_idx) = n + 1;
                        end
                    end
                end
            end
        end
    end
    
    x = Theta(:);
    y = Axis3(:);
    z = SNR(:);
    c = combined_metric(:);
    
    valid_idx = ~isnan(c);
    
    total_points = numel(combined_metric);
    non_overlapping_points = sum(valid_idx);
    
    fprintf('3D Phase Diagram Statistics:\n');
    fprintf('- Total parameter combinations: %d\n', total_points);
    fprintf('- Populated datapoints: %d\n', non_overlapping_points);
    fprintf('- Coverage: %.1f%%\n', 100 * non_overlapping_points / total_points);
    
    overlapping_idx = point_counts(:) > 1 & valid_idx;
    non_overlapping_idx = point_counts(:) == 1 & valid_idx;
    
    hold on;
    h = [];
    
    if any(non_overlapping_idx)
        h = scatter3(x(non_overlapping_idx), y(non_overlapping_idx), z(non_overlapping_idx), ...
            100, combined_metric(non_overlapping_idx), 'o', 'filled');
    end
    
    if any(overlapping_idx)
        h2 = scatter3(x(overlapping_idx), y(overlapping_idx), z(overlapping_idx), ...
            100, combined_metric(overlapping_idx), '^', 'filled');
        if isempty(h), h = h2; end
    end
    
    hold off;
    
    try
        cmap_data = load('metric_colormapv2.mat');
        colormap(cmap_data.CustomColormap);
    catch
        colormap('parula');
    end
    
    colorbar;
    clim([0 1]);
    
    xlabel('defect density');
    ylabel(axis3_label);
    zlabel('SNR');
    title(sprintf('%s (Combined %d Runs)', title_str, length(dataset_metrics_array)));
    
    grid on;
    view(-35, 79);
    
    cb = colorbar;
    ylabel(cb, 'Performance Score');
    
    ax = gca;
    ax.XScale = 'log';
    ax.YScale = 'linear';
    ax.ZScale = 'linear';
    
    set(ax, 'XDir', 'normal');
    set(ax, 'YDir', 'reverse');
    set(ax, 'ZDir', 'reverse');
    
    ax.XTick = unique_theta;
    ax.YTick = unique_axis3;
    
    if num_snr_values > 1
        snr_spaced = linspace(min(unique_snr), max(unique_snr), num_snr_values);
        ax.ZTick = snr_spaced;
        ax.ZTickLabel = arrayfun(@(x) sprintf('%.1f', x), unique_snr, 'UniformOutput', false);
    else
        ax.ZTick = unique_snr;
        ax.ZTickLabel = arrayfun(@(x) sprintf('%.1f', x), unique_snr, 'UniformOutput', false);
    end
    
    ax.XTickLabel = arrayfun(@(x) sprintf('%.0e', x), unique_theta, 'UniformOutput', false);
    ax.YTickLabel = arrayfun(@(x) sprintf(axis3_tick_fmt, x), unique_axis3, 'UniformOutput', false);
    
    ax.XTickLabelRotation = 90;
    axis square;
    ax.TickLength = [0.02 0.02];
    ax.FontSize = 15;
end

function avg_data = average_over_reps(data)
    if ndims(data) == 4
        avg_data = mean(data, 4, 'omitnan');
    else
        avg_data = data;
    end
end