function plot_defects_vs_kernel_similarity_by_snr(combined_metrics)
% For each unique SNR, plot a 2D scatter: x = avg number of defects, y = kernel similarity

SNRs = combined_metrics.SNR_values;
num_SNR = length(SNRs);

% Prepare data
num_points = numel(combined_metrics.kernel_quality_final);
defects = nan(num_points, 1);
kernel_sim = nan(num_points, 1);
SNR_idx = nan(num_points, 1);
defect_density_idx = nan(num_points, 1);

% The combined_metrics.kernel_quality_final array can be 3D [SNR, theta_cap, Nobs] 
% or 4D [SNR, theta_cap, Nobs, repetition]
% In both cases: idx1 = SNR index, idx2 = theta_cap (defect density) index, idx3 = Nobs index
% For 4D: idx4 = repetition index (we'll iterate over all repetitions)
array_dims = size(combined_metrics.kernel_quality_final);
num_dims = length(array_dims);
if num_dims < 3 || num_dims > 4
    error('Expected 3D or 4D array for kernel_quality_final: [SNR, theta_cap, Nobs] or [SNR, theta_cap, Nobs, repetition]');
end
if length(SNRs) ~= array_dims(1)
    warning('SNR_values length (%d) does not match array dimension 1 (%d)', length(SNRs), array_dims(1));
end

count = 0;
for i = 1:numel(combined_metrics.kernel_quality_final)
    kq = combined_metrics.kernel_quality_final(i);
    if isnan(kq)
        continue;
    end
    
    % Handle both 3D and 4D arrays
    if num_dims == 3
        [idx1, idx2, idx3] = ind2sub(array_dims, i);
        % idx1 = SNR index, idx2 = theta_cap (defect density) index, idx3 = Nobs index
        X0 = combined_metrics.X0{idx1, idx2, idx3};
    else % 4D
        [idx1, idx2, idx3, idx4] = ind2sub(array_dims, i);
        % idx1 = SNR index, idx2 = theta_cap (defect density) index, idx3 = Nobs index, idx4 = repetition index
        X0 = combined_metrics.X0{idx1, idx2, idx3, idx4};
    end
    
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
    defect_density_idx(count) = idx2;  % idx2 is the defect density (theta_cap) index
end

defects = defects(1:count);
kernel_sim = kernel_sim(1:count);
SNR_idx = SNR_idx(1:count);
defect_density_idx = defect_density_idx(1:count);

% Only keep points with kernel similarity > 0.84
valid = kernel_sim > 0.84;
defects = defects(valid);
kernel_sim = kernel_sim(valid);
SNR_idx = SNR_idx(valid);
defect_density_idx = defect_density_idx(valid);

% Get defect density values
if isfield(combined_metrics, 'theta_cap_values')
    defect_density_values = combined_metrics.theta_cap_values;
elseif isfield(combined_metrics, 'defect_density_values')
    defect_density_values = combined_metrics.defect_density_values;
else
    error(['Could not find defect density values in combined_metrics. ' ...
           'The combine_metrics_for_plotting function should store theta_cap_values. ' ...
           'Array dimension 2 is %d, which should match the number of unique theta_cap values.'], ...
           array_dims(2));
end

% Verify the length matches the array dimension
if length(defect_density_values) ~= array_dims(2)
    warning('theta_cap_values length (%d) does not match array dimension 2 (%d)', ...
            length(defect_density_values), array_dims(2));
end

% Print theta_cap values used in this dataset
fprintf('\n=== Theta_cap (Defect Density) Values in Dataset ===\n');
fprintf('Number of unique theta_cap values: %d\n', length(defect_density_values));
fprintf('Values:\n');
for i = 1:length(defect_density_values)
    fprintf('  [%d] %.6e\n', i, defect_density_values(i));
end
fprintf('====================================================\n\n');

% After preparing defects, kernel_sim, SNR_idx, etc.
x_min = min(defects);
x_max = max(defects);
% Add padding to x-axis limits (multiplicative for log scale)
padding_factor = 1.1;  % 10% padding on each side
x_min_padded = x_min / padding_factor;
x_max_padded = x_max * padding_factor;

% Map defect density indices to actual values
defect_density_actual = defect_density_values(defect_density_idx);

% Generate colors for each unique defect density value
num_densities = length(defect_density_values);
% Use a colormap to generate distinct colors for each density
cmap = parula(num_densities);
density_colors_map = containers.Map('KeyType', 'double', 'ValueType', 'any');
for i = 1:num_densities
    density_colors_map(defect_density_values(i)) = cmap(i, :);
end

% Plot 1: horizontal subplots, sharing the same Y axis
figure;
t = tiledlayout(1, num_SNR, 'TileSpacing', 'tight', 'Padding', 'compact');
axes_handles = gobjects(num_SNR, 1);  % Store axes handles

for s = 1:num_SNR
    ax = nexttile;
    axes_handles(s) = ax;
    
    % Plot points grouped by defect density with distinct colors
    hold(ax, 'on');
    legend_handles = [];
    legend_labels = {};
    
    for d = 1:num_densities
        density_val = defect_density_values(d);
        % Match by index to avoid floating point precision issues
        density_mask = defect_density_idx == d;
        snr_mask = SNR_idx == s;
        combined_mask = density_mask & snr_mask;
        
        if any(combined_mask)
            h = scatter(ax, defects(combined_mask), kernel_sim(combined_mask), 60, ...
                       density_colors_map(density_val), 'filled');
            legend_handles(end+1) = h;
            legend_labels{end+1} = sprintf('%.2e', density_val);
        end
    end
    hold(ax, 'off');
    
    xlabel('Average occurrence of Defects');
    xscale('log');
    % Set ylabel only on first subplot to share y-axis
    if s == 1
        ylabel('Kernel Similarity');
    else
        ylabel('');
        ax.YTickLabel = [];
    end
    title(sprintf('SNR = %.2f', SNRs(s)));
    grid on;
    % Set same limits for all subplots to ensure equal x-axis length
    ylim([0.84 1]);
    xlim([x_min_padded, x_max_padded]);
    
    % Show legend only on first subplot with font size 12 and title
    if s == 1 && ~isempty(legend_handles)
        lg = legend(ax, legend_handles, legend_labels, 'Location', 'best');
        lg.FontSize = 12;
        % Add title as text annotation since Title property may be read-only
        try
            lg.Title.String = 'Defect Density';
            lg.Title.FontSize = 12;
        catch
            % If Title is read-only, use text annotation instead
            lg_pos = lg.Position;
            text(ax, lg_pos(1) + lg_pos(3)/2, lg_pos(2) + lg_pos(4), ...
                 'Defect Density', 'HorizontalAlignment', 'center', ...
                 'FontSize', 12, 'FontWeight', 'bold', 'Units', 'normalized');
        end
    end
end

% Link axes to share y-axis and ensure same x-axis range
linkaxes(axes_handles, 'xy');  % Link both x and y axes

sgtitle('Defects vs Kernel Similarity by SNR (Colored by Defect Density)');
end 