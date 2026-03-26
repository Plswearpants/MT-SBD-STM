function metrics2heat_properGen(dataset_metrics, mode)
    % Create 3D heatspace visualization for both kernel similarity and combined activation metrics
    % Input:
    %   dataset_metrics: structure containing:
    %       - SNR_values
    %       - theta_cap_values
    %       - Nobs_values (mode 1)
    %       - side_length_ratio_values (mode 2, optional)
    %       - repetition_values (optional, for 4D arrays)
    %       - kernel_similarity [SNR × theta_cap × Nobs] or [SNR × theta_cap × Nobs × rep]
    %       - activation_metrics [SNR × theta_cap × Nobs] or [SNR × theta_cap × Nobs × rep]
    % Note: If metrics have a 4th dimension (repetitions), they will be averaged over repetitions
    %   mode: optional plotting mode
    %       - "1" (default): [SNR × theta_cap × Nobs]
    %       - "2": [SNR × theta_cap × side_length_ratio]

    if nargin < 2 || isempty(mode)
        mode = "1";
    else
        mode = string(mode);
    end

    switch mode
        case "1"
            axis3_values = dataset_metrics.Nobs_values;
            axis3_label = 'N_{obs}';
            axis3_tick_fmt = '%.0f';
            kernel_data = dataset_metrics.kernel_quality_final;
            combined_data = dataset_metrics.combined_activationScore;
        case "2"
            if isfield(dataset_metrics, 'side_length_ratio_values')
                axis3_values = dataset_metrics.side_length_ratio_values;
            elseif isfield(dataset_metrics, 'area_ratio_values')
                % Backward-compatible fallback if older field name is used
                axis3_values = dataset_metrics.area_ratio_values;
            else
                error(['mode "2" requires dataset_metrics.side_length_ratio_values ' ...
                    '(or area_ratio_values for backward compatibility).']);
            end
            axis3_label = 'side length ratio';
            axis3_tick_fmt = '%.3f';
            if isfield(dataset_metrics, 'kernel_quality_final_by_side_length_ratio') && ...
               isfield(dataset_metrics, 'combined_activationScore_by_side_length_ratio')
                kernel_data = dataset_metrics.kernel_quality_final_by_side_length_ratio;
                combined_data = dataset_metrics.combined_activationScore_by_side_length_ratio;
            else
                % Fallback for older loaders: interpret existing tensors with mode-2 axis.
                kernel_data = dataset_metrics.kernel_quality_final;
                combined_data = dataset_metrics.combined_activationScore;
            end
        otherwise
            error('Unsupported mode "%s". Use "1" or "2".', mode);
    end
    
    % Average over repetitions if 4D arrays are present
    kernel_quality_avg = average_over_repetitions(kernel_data);
    combined_activation_avg = average_over_repetitions(combined_data);
    
    % Create figure with two subplots
    figure('Position', [100 100 1600 800]);
    
    % 1. Kernel Similarity Heatspace
    subplot(1,2,1);
    h1 = plot_3D_heatspace(dataset_metrics.SNR_values, ...
                          dataset_metrics.theta_cap_values, ...
                          axis3_values, ...
                          kernel_quality_avg, ...
                          'Kernel Similarity (Averaged over Repetitions)', ...
                          axis3_label, axis3_tick_fmt);
    
    % 2. Combined Activation Score Heatspace
    subplot(1,2,2);
    h2 = plot_3D_heatspace(dataset_metrics.SNR_values, ...
                          dataset_metrics.theta_cap_values, ...
                          axis3_values, ...
                          combined_activation_avg, ...
                          'Combined Activation Score (Averaged over Repetitions)', ...
                          axis3_label, axis3_tick_fmt);
end

function h = plot_3D_heatspace(SNR_values, theta_cap_values, axis3_values, metric_values, title_str, axis3_label, axis3_tick_fmt)
    % Create meshgrid matching the data structure: Metric(SNR, theta_cap, axis3)
    % The data structure is [SNR, theta_cap, axis3], so meshgrid should match this
    [SNR, Theta, Axis3] = meshgrid(SNR_values, theta_cap_values, axis3_values);
    
    % Permute the meshgrid to match the data structure dimensions [SNR, theta_cap, Nobs]
    % Original meshgrid creates [theta_cap, SNR, axis3], so we need to permute
    SNR = permute(SNR, [2, 1, 3]);  % Move SNR to first dimension
    Theta = permute(Theta, [2, 1, 3]);  % Move theta_cap to second dimension
    Axis3 = permute(Axis3, [2, 1, 3]);  % Keep axis3 in third dimension
    
    % Remap SNR coordinates to evenly spaced positions
    num_snr_values = length(SNR_values);
    if num_snr_values > 1
        % Create evenly spaced SNR positions
        snr_spaced = linspace(min(SNR_values), max(SNR_values), num_snr_values);
        
        % Create mapping from actual SNR values to evenly spaced positions
        snr_mapping = containers.Map(SNR_values, snr_spaced);
        
        % Remap SNR coordinates
        SNR_remapped = zeros(size(SNR));
        for i = 1:numel(SNR)
            SNR_remapped(i) = snr_mapping(SNR(i));
        end
        SNR = SNR_remapped;
    end
    
    % Ensure all inputs are column vectors
    x = Theta(:);   % Column vector for x coordinates (theta_cap)
    y = Axis3(:);   % Column vector for y coordinates (axis3)
    z = SNR(:);     % Column vector for z coordinates (SNR) - now evenly spaced
    c = metric_values(:);  % Column vector for colors

    % Create scatter plot with color based on metric values
    h = scatter3(x, y, z, 100, c, 'filled');
    
    % Load custom colormap
    cmap_data = load('metric_colormapv3.mat');
    colormap(cmap_data.CustomColormap);
    
    % Customize appearance
    colorbar;
    clim([0 1]);
    
    % Add labels and title
    xlabel('defect density');
    ylabel(axis3_label);
    zlabel('SNR');
    title(title_str);
    
    % Add grid and adjust view
    grid on;
    view(-35, 85);
    
    % Add colorbar label
    c = colorbar;
    ylabel(c, 'Performance Score');
    
    % Make axis labels more readable and set directions
    ax = gca;
    ax.XScale = 'log';  % Theta cap in log scale
    ax.YScale = 'linear';  % axis3 in linear scale
    ax.ZScale = 'linear';  % SNR in linear scale
    
    % Set axis directions
    set(ax, 'XDir', 'normal');   % theta_cap: low -> high
    set(ax, 'YDir', 'reverse');  % axis3: high -> low
    set(ax, 'ZDir', 'reverse');  % SNR: high -> low
    
    % Set tick values to match the actual data points
    ax.XTick = theta_cap_values;
    ax.YTick = axis3_values;
    
    % Set SNR tick positions to the evenly spaced positions
    if num_snr_values > 1
        snr_spaced = linspace(min(SNR_values), max(SNR_values), num_snr_values);
        ax.ZTick = snr_spaced;
    else
        ax.ZTick = SNR_values;
    end
    
    % Format tick labels
    ax.XTickLabel = arrayfun(@(x) sprintf('%.0e', x), theta_cap_values, 'UniformOutput', false);
    ax.YTickLabel = arrayfun(@(x) sprintf(axis3_tick_fmt, x), axis3_values, 'UniformOutput', false);
    
    % Format SNR tick labels to show actual SNR values at evenly spaced positions
    if num_snr_values > 1
        % Create tick labels that correspond to the actual SNR values
        snr_tick_labels = arrayfun(@(x) sprintf('%.1f', x), SNR_values, 'UniformOutput', false);
        ax.ZTickLabel = snr_tick_labels;
    else
        ax.ZTickLabel = arrayfun(@(x) sprintf('%.1f', x), SNR_values, 'UniformOutput', false);
    end
    
    % Rotate x-axis labels for better readability
    ax.XTickLabelRotation = 45;
end

function avg_data = average_over_repetitions(data)
    % Average over the 4th dimension (repetitions) if present
    % Input: data can be 3D [SNR × theta × N_obs] or 4D [SNR × theta × N_obs × rep]
    % Output: 3D array [SNR × theta × N_obs] averaged over repetitions
    
    if ndims(data) == 4
        % Average over the 4th dimension, ignoring NaN values
        avg_data = mean(data, 4, 'omitnan');
    elseif ndims(data) == 3
        % Already 3D, return as is
        avg_data = data;
    else
        % Unexpected dimensions, return as is
        avg_data = data;
    end
end
