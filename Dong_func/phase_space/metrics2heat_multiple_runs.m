function metrics2heat_multiple_runs(dataset_metrics_array, run_names)
    % Create 3D heatspace visualization for multiple experiment runs
    % Inputs:
    %   dataset_metrics_array: Cell array of dataset_metrics structures
    %   run_names: Cell array of names for each run (optional)
    
    if nargin < 2
        run_names = cell(length(dataset_metrics_array), 1);
        for i = 1:length(dataset_metrics_array)
            run_names{i} = sprintf('Run %d', i);
        end
    end
    
    % Validate inputs
    if length(dataset_metrics_array) ~= length(run_names)
        error('Number of dataset_metrics must match number of run_names');
    end
    
    % Create figure with two subplots
    figure('Position', [100 100 1600 800]);
    
    % 1. Kernel Similarity Heatspace
    subplot(1,2,1);
    h1 = plot_3D_heatspace_multiple(dataset_metrics_array, run_names, 'kernel_quality_final', 'Kernel Similarity');
    
    % 2. Combined Activation Score Heatspace
    subplot(1,2,2);
    h2 = plot_3D_heatspace_multiple(dataset_metrics_array, run_names, 'combined_activationScore', 'Combined Activation Score');
end

function h = plot_3D_heatspace_multiple(dataset_metrics_array, run_names, metric_field, title_str)
    % Plot multiple runs with averaging of overlapping points
    
    % Collect all unique parameter values across all runs
    all_snr_values = [];
    all_theta_values = [];
    all_nobs_values = [];
    
    for i = 1:length(dataset_metrics_array)
        metrics = dataset_metrics_array{i};
        all_snr_values = [all_snr_values; metrics.SNR_values(:)];
        all_theta_values = [all_theta_values; metrics.theta_cap_values(:)];
        all_nobs_values = [all_nobs_values; metrics.Nobs_values(:)];
    end
    
    % Get unique sorted values
    unique_snr = unique(all_snr_values);
    unique_theta = unique(all_theta_values);
    unique_nobs = unique(all_nobs_values);
    
    % Create combined parameter space
    [SNR, Theta, Nobs] = meshgrid(unique_snr, unique_theta, unique_nobs);
    
    % Permute to match data structure dimensions [SNR, theta_cap, Nobs]
    SNR = permute(SNR, [2, 1, 3]);
    Theta = permute(Theta, [2, 1, 3]);
    Nobs = permute(Nobs, [2, 1, 3]);
    
    % Remap SNR coordinates to evenly spaced positions
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
    
    % Initialize combined metric array
    combined_metric = nan(size(SNR));
    point_counts = zeros(size(SNR));  % Track how many runs contribute to each point
    
    % Combine metrics from all runs
    for run_idx = 1:length(dataset_metrics_array)
        metrics = dataset_metrics_array{run_idx};
        
        % Get metric values for this run
        if isfield(metrics, metric_field)
            run_metric = metrics.(metric_field);
        else
            warning('Field %s not found in run %d', metric_field, run_idx);
            continue;
        end
        
        % Map run's parameter space to combined parameter space
        for i = 1:length(metrics.SNR_values)
            for j = 1:length(metrics.theta_cap_values)
                for k = 1:length(metrics.Nobs_values)
                    % Find corresponding indices in combined space
                    [~, snr_idx] = min(abs(unique_snr - metrics.SNR_values(i)));
                    [~, theta_idx] = min(abs(unique_theta - metrics.theta_cap_values(j)));
                    [~, nobs_idx] = min(abs(unique_nobs - metrics.Nobs_values(k)));
                    
                    % Check if this point has a valid metric value
                    if ~isnan(run_metric(i, j, k))
                        if isnan(combined_metric(snr_idx, theta_idx, nobs_idx))
                            combined_metric(snr_idx, theta_idx, nobs_idx) = run_metric(i, j, k);
                            point_counts(snr_idx, theta_idx, nobs_idx) = 1;
                        else
                            % Average with existing value
                            total = combined_metric(snr_idx, theta_idx, nobs_idx) * point_counts(snr_idx, theta_idx, nobs_idx);
                            total = total + run_metric(i, j, k);
                            point_counts(snr_idx, theta_idx, nobs_idx) = point_counts(snr_idx, theta_idx, nobs_idx) + 1;
                            combined_metric(snr_idx, theta_idx, nobs_idx) = total / point_counts(snr_idx, theta_idx, nobs_idx);
                        end
                    end
                end
            end
        end
    end
    
    % Prepare data for plotting
    x = Theta(:);
    y = Nobs(:);
    z = SNR(:);
    c = combined_metric(:);
    
    % Remove NaN values
    valid_idx = ~isnan(c);
    x = x(valid_idx);
    y = y(valid_idx);
    z = z(valid_idx);
    c = c(valid_idx);
    
    % Print statistics about non-overlapping datapoints
    total_points = numel(combined_metric);
    non_overlapping_points = sum(valid_idx);
    overlapping_points = total_points - non_overlapping_points;
    
    fprintf('3D Phase Diagram Statistics:\n');
    fprintf('- Total parameter combinations: %d\n', total_points);
    fprintf('- Non-overlapping datapoints: %d\n', non_overlapping_points);
    fprintf('- Overlapping datapoints (averaged): %d\n', overlapping_points);
    fprintf('- Coverage: %.1f%%\n', 100 * non_overlapping_points / total_points);
    
    % Separate overlapping and non-overlapping points
    overlapping_idx = point_counts(:) > 1;
    overlapping_idx = overlapping_idx & valid_idx;
    non_overlapping_idx = point_counts(:) == 1;
    non_overlapping_idx = non_overlapping_idx & valid_idx;
    
    % Create scatter plots with different markers
    hold on;
    
    % Plot non-overlapping points with circles
    if any(non_overlapping_idx)
        x_non = Theta(:);
        y_non = Nobs(:);
        z_non = SNR(:);
        c_non = combined_metric(:);
        
        x_non = x_non(non_overlapping_idx);
        y_non = y_non(non_overlapping_idx);
        z_non = z_non(non_overlapping_idx);
        c_non = c_non(non_overlapping_idx);
        
        h1 = scatter3(x_non, y_non, z_non, 100, c_non, 'o', 'filled');
    end
    
    % Plot overlapping points with triangles
    if any(overlapping_idx)
        x_over = Theta(:);
        y_over = Nobs(:);
        z_over = SNR(:);
        c_over = combined_metric(:);
        
        x_over = x_over(overlapping_idx);
        y_over = y_over(overlapping_idx);
        z_over = z_over(overlapping_idx);
        c_over = c_over(overlapping_idx);
        
        h2 = scatter3(x_over, y_over, z_over, 100, c_over, '^', 'filled');
    end
    
    hold off;
    
    % Return the first handle for compatibility
    if any(non_overlapping_idx)
        h = h1;
    elseif any(overlapping_idx)
        h = h2;
    else
        h = [];
    end
    
    % Load custom colormap
    cmap_data = load('metric_colormapv2.mat');
    colormap(cmap_data.CustomColormap);
    
    % Customize appearance
    colorbar;
    clim([0 1]);
    
    % Add labels and title
    xlabel('');
    ylabel('');
    zlabel('');
    title(sprintf('%s (Combined %d Runs)', title_str, length(dataset_metrics_array)));
    
    % Add grid and adjust view
    grid on;
    view(-35, 79);
    
    % Add colorbar label
    c = colorbar;
    ylabel(c, 'Performance Score');
    
    % Make axis labels more readable and set directions
    ax = gca;
    ax.XScale = 'log';
    ax.YScale = 'linear';
    ax.ZScale = 'linear';
    
    % Set axis directions
    set(ax, 'XDir', 'normal');
    set(ax, 'YDir', 'reverse');
    set(ax, 'ZDir', 'reverse');
    
    % Set tick values
    ax.XTick = unique_theta;
    ax.YTick = unique_nobs;
    
    if num_snr_values > 1
        snr_spaced = linspace(min(unique_snr), max(unique_snr), num_snr_values);
        ax.ZTick = snr_spaced;
    else
        ax.ZTick = unique_snr;
    end
    
    % Format tick labels
    ax.XTickLabel = arrayfun(@(x) sprintf('%.0e', x), unique_theta, 'UniformOutput', false);
    ax.YTickLabel = arrayfun(@(x) sprintf('%.0f', x), unique_nobs, 'UniformOutput', false);
    
    if num_snr_values > 1
        snr_tick_labels = arrayfun(@(x) sprintf('%.1f', x), unique_snr, 'UniformOutput', false);
        ax.ZTickLabel = snr_tick_labels;
    else
        ax.ZTickLabel = arrayfun(@(x) sprintf('%.1f', x), unique_snr, 'UniformOutput', false);
    end
    
    % Rotate x-axis labels
    ax.XTickLabelRotation = 90;
    axis square;
    
    % Increase tick mark size
    ax.TickLength = [0.02 0.02];  % [length width] - default is [0.01 0.025]
    
    % Increase tick label font size
    ax.FontSize = 15;  % Increase from default (usually 10-12)
end 