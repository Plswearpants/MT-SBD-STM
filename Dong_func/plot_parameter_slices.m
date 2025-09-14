function plot_parameter_slices(metrics)
    % Define metrics to plot (replace activation and demixing scores with combined score)
    metric_fields = {'kernel_quality_final', 'combined_score', 'runtime'};
    metric_names = {'Kernel Quality', 'Combined Score', 'Runtime (s)'};

    % Compute combined score if not already present
    if ~isfield(metrics, 'combined_score')
        metrics.combined_score = computeCombined_activationScore(metrics.demixing_score, ...
                                                metrics.activation_accuracy_final);
    end
    
    % Get parameter values
    lambda1_values = metrics.lambda1_values;
    mini_loop_values = metrics.mini_loop_values;
    num_datasets = metrics.num_datasets;
    
    % Create two figures: one for each parameter
    % 1. Lambda1 variation (fixed mini_loop)
    fig1 = figure('Position', [100 100 1200 800], ...
        'Name', 'Parameter Slice: Lambda1 Variation');
    
    % Add controls for lambda1 figure
    control_panel1 = uipanel('Parent', fig1, ...
        'Units', 'normalized', ...
        'Position', [0 0.95 1 0.05], ...
        'BackgroundColor', get(fig1, 'Color'));
    
    % Mini-loop selector for lambda1 plot
    uicontrol('Parent', control_panel1, ...
        'Style', 'text', ...
        'Units', 'normalized', ...
        'Position', [0.02 0.1 0.1 0.8], ...
        'String', 'Mini-loop value:', ...
        'BackgroundColor', get(fig1, 'Color'));
    
    mini_loop_selector = uicontrol('Parent', control_panel1, ...
        'Style', 'popup', ...
        'Units', 'normalized', ...
        'Position', [0.13 0.1 0.1 0.8], ...
        'String', arrayfun(@(x) sprintf('Mini-loop=%d', x), ...
        mini_loop_values, 'UniformOutput', false), ...
        'Callback', @updateLambda1Plot);
    
    % Dataset highlight controls
    addHighlightControls(control_panel1, metrics, @updateLambda1Plot);
    
    % Create plot panel
    plot_panel1 = uipanel('Parent', fig1, ...
        'Units', 'normalized', ...
        'Position', [0 0 1 0.95]);
    
    % 2. Mini-loop variation (fixed lambda1)
    fig2 = figure('Position', [150 150 1200 800], ...
        'Name', 'Parameter Slice: Mini-loop Variation');
    
    % Add controls for mini-loop figure
    control_panel2 = uipanel('Parent', fig2, ...
        'Units', 'normalized', ...
        'Position', [0 0.95 1 0.05], ...
        'BackgroundColor', get(fig2, 'Color'));
    
    % Lambda1 selector for mini-loop plot
    uicontrol('Parent', control_panel2, ...
        'Style', 'text', ...
        'Units', 'normalized', ...
        'Position', [0.02 0.1 0.1 0.8], ...
        'String', 'Lambda1 value:', ...
        'BackgroundColor', get(fig2, 'Color'));
    
    lambda1_selector = uicontrol('Parent', control_panel2, ...
        'Style', 'popup', ...
        'Units', 'normalized', ...
        'Position', [0.13 0.1 0.1 0.8], ...
        'String', arrayfun(@(x) sprintf('Lambda1=%.3f', x), ...
        lambda1_values, 'UniformOutput', false), ...
        'Callback', @updateMiniLoopPlot);
    
    % Dataset highlight controls
    addHighlightControls(control_panel2, metrics, @updateMiniLoopPlot);
    
    % Create plot panel
    plot_panel2 = uipanel('Parent', fig2, ...
        'Units', 'normalized', ...
        'Position', [0 0 1 0.95]);
    
    % Store data in figures
    fig1.UserData.metrics = metrics;
    fig1.UserData.mini_loop_selector = mini_loop_selector;
    fig1.UserData.highlighted_dataset = 1;
    fig1.UserData.other_alpha = 0.3;
    fig1.UserData.plot_panel = plot_panel1;
    
    fig2.UserData.metrics = metrics;
    fig2.UserData.lambda1_selector = lambda1_selector;
    fig2.UserData.highlighted_dataset = 1;
    fig2.UserData.other_alpha = 0.3;
    fig2.UserData.plot_panel = plot_panel2;
    
    % Initial plots
    updateLambda1Plot();
    updateMiniLoopPlot();
    
    function updateLambda1Plot(~,~)
        mini_loop_idx = fig1.UserData.mini_loop_selector.Value;
        highlighted_dataset = fig1.UserData.highlighted_dataset;
        other_alpha = fig1.UserData.other_alpha;
        
        % Clear previous plots
        delete(fig1.UserData.plot_panel.Children);
        
        % Create new tiledlayout with 1x3 format
        t = tiledlayout(fig1.UserData.plot_panel, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
        
        % Plot each metric
        for m = 1:length(metric_fields)
            ax = nexttile(t);
            plotMetricSlice(ax, metrics, metric_fields{m}, metric_names{m}, ...
                'lambda1', mini_loop_idx, highlighted_dataset, other_alpha);
        end
        
        sgtitle(sprintf('Metrics vs Lambda1 (Mini-loop=%d)', mini_loop_values(mini_loop_idx)));
    end
    
    function updateMiniLoopPlot(~,~)
        lambda1_idx = fig2.UserData.lambda1_selector.Value;
        highlighted_dataset = fig2.UserData.highlighted_dataset;
        other_alpha = fig2.UserData.other_alpha;
        
        % Clear previous plots
        delete(fig2.UserData.plot_panel.Children);
        
        % Create new tiledlayout with 1x3 format
        t = tiledlayout(fig2.UserData.plot_panel, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
        
        % Plot each metric
        for m = 1:length(metric_fields)
            ax = nexttile(t);
            plotMetricSlice(ax, metrics, metric_fields{m}, metric_names{m}, ...
                'mini_loop', lambda1_idx, highlighted_dataset, other_alpha);
        end
        
        sgtitle(sprintf('Metrics vs Mini-loop (Lambda1=%.3f)', lambda1_values(lambda1_idx)));
    end
end

function addHighlightControls(panel, metrics, callback_fn)
    % Add dataset selector
    uicontrol('Parent', panel, ...
        'Style', 'text', ...
        'Units', 'normalized', ...
        'Position', [0.24 0.1 0.1 0.8], ...
        'String', 'Highlight Dataset:', ...
        'BackgroundColor', get(panel.Parent, 'Color'));
    
    uicontrol('Parent', panel, ...
        'Style', 'popup', ...
        'Units', 'normalized', ...
        'Position', [0.35 0.1 0.1 0.8], ...
        'String', arrayfun(@(x) sprintf('Dataset %d', x), ...
        1:metrics.num_datasets, 'UniformOutput', false), ...
        'Value', 1, ...
        'Callback', @(~,~) updateHighlight(panel, callback_fn));
    
    % Add transparency slider
    uicontrol('Parent', panel, ...
        'Style', 'text', ...
        'Units', 'normalized', ...
        'Position', [0.46 0.1 0.12 0.8], ...
        'String', 'Other Datasets Alpha:', ...
        'BackgroundColor', get(panel.Parent, 'Color'));
    
    uicontrol('Parent', panel, ...
        'Style', 'slider', ...
        'Units', 'normalized', ...
        'Position', [0.59 0.1 0.1 0.8], ...
        'Min', 0, 'Max', 1, 'Value', 0.3, ...
        'Callback', @(~,~) updateHighlight(panel, callback_fn));
    
    % Add description text box
    if isfield(metrics, 'dataset_descriptions')
        uicontrol('Parent', panel, ...
            'Style', 'text', ...
            'Units', 'normalized', ...
            'Position', [0.7 0.1 0.28 0.8], ...
            'String', metrics.dataset_descriptions{1}, ...
            'BackgroundColor', get(panel.Parent, 'Color'), ...
            'HorizontalAlignment', 'left');
    end
end

function updateHighlight(panel, callback_fn)
    fig = panel.Parent;
    controls = panel.Children;
    
    % Find dataset selector and slider
    dataset_selector = findobj(controls, 'Style', 'popup');
    alpha_slider = findobj(controls, 'Style', 'slider');
    description_text = findobj(controls, 'Style', 'text', 'Position', [0.7 0.1 0.28 0.8]);
    
    % Update stored values
    fig.UserData.highlighted_dataset = dataset_selector.Value;
    fig.UserData.other_alpha = alpha_slider.Value;
    
    % Update description if available
    try
        % Get the current dataset index
        current_idx = dataset_selector.Value;
        
        % Check if we have valid descriptions
        if ~isempty(description_text) && ...
           isfield(fig.UserData.metrics, 'dataset_descriptions') && ...
           iscell(fig.UserData.metrics.dataset_descriptions) && ...
           current_idx <= length(fig.UserData.metrics.dataset_descriptions)
            
            % Get the description for current dataset
            current_description = fig.UserData.metrics.dataset_descriptions{current_idx};
            
            % Update text if description exists
            if ~isempty(current_description)
                description_text.String = current_description;
            else
                description_text.String = sprintf('Dataset %d', current_idx);
            end
        else
            description_text.String = sprintf('Dataset %d', current_idx);
        end
    catch
        % Fallback if anything goes wrong
        if ~isempty(description_text)
            description_text.String = sprintf('Dataset %d', dataset_selector.Value);
        end
    end
    
    % Update plots
    callback_fn();
end

function plotMetricSlice(ax, metrics, metric_field, metric_name, vary_param, fixed_idx, highlighted_dataset, other_alpha)
    data = metrics.(metric_field);  % 3D array: datasets × mini_loop × lambda1
    
    if strcmp(vary_param, 'lambda1')
        x_values = metrics.lambda1_values;
        x_label = '\lambda_1';
        y_values = squeeze(data(:, fixed_idx, :));  % Get all lambda1 values for fixed mini_loop
    else
        x_values = metrics.mini_loop_values;
        x_label = 'mini\_loop';
        y_values = squeeze(data(:, :, fixed_idx));  % Get all mini_loop values for fixed lambda1
    end
    
    hold(ax, 'on');
    
    % Plot each dataset
    for i = 1:metrics.num_datasets
        % Set transparency based on highlighting
        if i == highlighted_dataset
            alpha = 1;
            line_width = 2;
        else
            alpha = other_alpha;
            line_width = 1;
        end
        
        plot(ax, x_values, y_values(i,:), 'o-', ...
            'LineWidth', line_width, ...
            'MarkerSize', 6, ...
            'Color', [get_dataset_color(i), alpha]);
    end
    
    xlabel(ax, x_label);
    ylabel(ax, metric_name);
    grid(ax, 'on');
    
    if strcmp(vary_param, 'lambda1')
        set(ax, 'XScale', 'log');
    end
    
    hold(ax, 'off');
end

function color = get_dataset_color(idx)
    colors = [
        0.8500    0.3250    0.0980;  % Orange
        0         0.4470    0.7410;  % Blue
        0.9290    0.6940    0.1250;  % Yellow
        0.4940    0.1840    0.5560;  % Purple
        0.4660    0.6740    0.1880;  % Green
        0.6350    0.0780    0.1840;  % Dark Red
        0         0.7500    0.7500;  % Cyan
        0.7500    0         0.7500;  % Magenta
        0.2500    0.2500    0.2500;  % Gray
        0.9500    0.5000    0.2000;  % Light Orange
        0.1000    0.5000    0.5000;  % Teal
        0.5000    0.5000    0;       % Olive
    ];
    color = colors(mod(idx-1, size(colors,1)) + 1, :);
end

function plotMetricOverview(ax, metrics, metric_field, metric_name, x_param)
    data = metrics.(metric_field);  % 3D array: datasets × mini_loop × lambda1
    markers = {'o', 's', 'd', '^', 'v', '>', '<', 'p', 'h'};
    line_alpha = 0.1;  % Transparency for connecting lines
    
    if strcmp(x_param, 'lambda1')
        x_values = metrics.lambda1_values;
        x_label = '\lambda_1';
        legend_values = metrics.mini_loop_values;
        legend_prefix = 'Mini-loop=';
        
        % Plot for each mini_loop value
        for j = 1:length(legend_values)
            marker_idx = mod(j-1, length(markers)) + 1;
            
            % Get data for all datasets at this mini_loop value
            y_values = squeeze(data(:, j, :));  % datasets × lambda1
            
            % First plot connecting lines with low alpha
            for i = 1:metrics.num_datasets
                plot(ax, x_values, y_values(i,:), '-', ...
                    'Color', [0 0 0, line_alpha], ...
                    'HandleVisibility', 'off');
            end
            
            % Then plot scatter points
            [X, Y] = meshgrid(x_values, 1:metrics.num_datasets);
            scatter(ax, X(:), y_values(:), 50, markers{marker_idx}, ...
                'DisplayName', sprintf('%s%g', legend_prefix, legend_values(j)));
        end
        
        set(ax, 'XScale', 'log');
    else
        x_values = metrics.mini_loop_values;
        x_label = 'mini\_loop';
        legend_values = metrics.lambda1_values;
        legend_prefix = '\lambda_1=';
        
        for j = 1:length(legend_values)
            marker_idx = mod(j-1, length(markers)) + 1;
            
            % Get data for all datasets at this lambda1 value
            y_values = squeeze(data(:, :, j));  % datasets × mini_loop
            
            % First plot connecting lines with low alpha
            for i = 1:metrics.num_datasets
                plot(ax, x_values, y_values(i,:), '-', ...
                    'Color', [0 0 0, line_alpha], ...
                    'HandleVisibility', 'off');
            end
            
            % Then plot scatter points
            [X, Y] = meshgrid(x_values, 1:metrics.num_datasets);
            scatter(ax, X(:), y_values(:), 50, markers{marker_idx}, ...
                'DisplayName', sprintf('%s%.3f', legend_prefix, legend_values(j)));
        end
    end
    
    xlabel(ax, x_label);
    ylabel(ax, metric_name);
    grid(ax, 'on');
    legend('show', 'Location', 'eastoutside');
    hold(ax, 'off');
end