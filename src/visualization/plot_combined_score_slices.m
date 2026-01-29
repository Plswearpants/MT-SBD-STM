function plot_combined_score_slices(metrics)
    % Create interactive visualizations of combined score vs parameters
    % Similar to plot_parameter_slices.m but focused on combined score
    
    % Get parameter values
    lambda1_values = metrics.lambda1_values;
    mini_loop_values = metrics.mini_loop_values;
    num_datasets = metrics.num_datasets;
    
    % Compute combined score for all parameter combinations
    combined_scores = computeCombined_activationScore(metrics.demixing_score, ...
                                         metrics.activation_accuracy_final);
    
    % 1. Lambda1 variation (fixed mini_loop)
    fig1 = figure('Position', [100 100 800 600], ...
        'Name', 'Combined Score: Lambda1 Variation');
    
    % Add controls for lambda1 figure
    control_panel1 = uipanel('Parent', fig1, ...
        'Units', 'normalized', ...
        'Position', [0 0.85 1 0.15], ...
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
        'Position', [0 0 1 0.85]);
    
    % 2. Mini-loop variation (fixed lambda1)
    fig2 = figure('Position', [150 150 800 600], ...
        'Name', 'Combined ActivationScore: Mini-loop Variation');
    
    % Add controls for mini-loop figure
    control_panel2 = uipanel('Parent', fig2, ...
        'Units', 'normalized', ...
        'Position', [0 0.85 1 0.15], ...
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
        'Position', [0 0 1 0.85]);
    
    % Store data in figures
    fig1.UserData.metrics = metrics;
    fig1.UserData.combined_scores = combined_scores;
    fig1.UserData.mini_loop_selector = mini_loop_selector;
    fig1.UserData.highlighted_dataset = 1;
    fig1.UserData.other_alpha = 0.3;
    fig1.UserData.plot_panel = plot_panel1;
    
    fig2.UserData.metrics = metrics;
    fig2.UserData.combined_scores = combined_scores;
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
        
        % Create axes
        ax = axes('Parent', fig1.UserData.plot_panel);
        
        % Get data for current mini_loop
        y_values = squeeze(combined_scores(:, mini_loop_idx, :));  % datasets × lambda1
        
        hold(ax, 'on');
        
        % Plot each dataset
        for i = 1:num_datasets
            % Set transparency based on highlighting
            if i == highlighted_dataset
                alpha = 1;
                line_width = 2;
            else
                alpha = other_alpha;
                line_width = 1;
            end
            
            plot(ax, lambda1_values, y_values(i,:), 'o-', ...
                'LineWidth', line_width, ...
                'MarkerSize', 6, ...
                'Color', [get_dataset_color(i), alpha]);
        end
        
        xlabel(ax, '\lambda_1');
        ylabel(ax, 'Combined Activation Score');
        title(sprintf('Combined Activation Score vs Lambda1 (Mini-loop=%d)', mini_loop_values(mini_loop_idx)));
        grid(ax, 'on');
        set(ax, 'XScale', 'log');
        ylim([0 1]);
        
        hold(ax, 'off');
    end
    
    function updateMiniLoopPlot(~,~)
        lambda1_idx = fig2.UserData.lambda1_selector.Value;
        highlighted_dataset = fig2.UserData.highlighted_dataset;
        other_alpha = fig2.UserData.other_alpha;
        
        % Clear previous plots
        delete(fig2.UserData.plot_panel.Children);
        
        % Create axes
        ax = axes('Parent', fig2.UserData.plot_panel);
        
        % Get data for current lambda1
        y_values = squeeze(combined_scores(:, :, lambda1_idx));  % datasets × mini_loop
        
        hold(ax, 'on');
        
        % Plot each dataset
        for i = 1:num_datasets
            % Set transparency based on highlighting
            if i == highlighted_dataset
                alpha = 1;
                line_width = 2;
            else
                alpha = other_alpha;
                line_width = 1;
            end
            
            plot(ax, mini_loop_values, y_values(i,:), 'o-', ...
                'LineWidth', line_width, ...
                'MarkerSize', 6, ...
                'Color', [get_dataset_color(i), alpha]);
        end
        
        xlabel(ax, 'mini\_loop');
        ylabel(ax, 'Combined Activation Score');
        title(sprintf('Combined Activation Score vs Mini-loop (Lambda1=%.3f)', lambda1_values(lambda1_idx)));
        grid(ax, 'on');
        ylim([0 1]);
        
        hold(ax, 'off');
    end
end

% Helper functions from original implementation
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
end

function updateHighlight(panel, callback_fn)
    fig = panel.Parent;
    controls = panel.Children;
    
    % Find dataset selector and slider
    dataset_selector = findobj(controls, 'Style', 'popup');
    alpha_slider = findobj(controls, 'Style', 'slider');
    
    % Update stored values
    fig.UserData.highlighted_dataset = dataset_selector.Value;
    fig.UserData.other_alpha = alpha_slider.Value;
    
    % Update plots
    callback_fn();
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