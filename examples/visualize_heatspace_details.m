function visualize_heatspace_details(dataset_metrics)
    % First show heatspace for reference
    metrics2heatspace(dataset_metrics);
    
    % Get user input for parameters
    fprintf('\nAvailable parameter ranges:\n');
    fprintf('Theta Cap: [%.2e to %.2e]\n', min(dataset_metrics.theta_cap_values), max(dataset_metrics.theta_cap_values));
    fprintf('Area Ratio: [%.3f to %.3f]\n', min(dataset_metrics.area_ratio_values), max(dataset_metrics.area_ratio_values));
    fprintf('SNR: [%.1f to %.1f]\n', min(dataset_metrics.SNR_values), max(dataset_metrics.SNR_values));
    
    % Get user input
    theta_cap = input('Enter Theta Cap value: ');
    area_ratio = input('Enter Area Ratio value: ');
    snr = input('Enter SNR value: ');
    
    % Find nearest indices in the parameter space
    [~, theta_idx] = min(abs(dataset_metrics.theta_cap_values - theta_cap));
    [~, area_idx] = min(abs(dataset_metrics.area_ratio_values - area_ratio));
    [~, snr_idx] = min(abs(dataset_metrics.SNR_values - snr));
    
    % Print actual values being used (in case of rounding)
    fprintf('\nUsing nearest available parameters:\n');
    fprintf('Theta Cap: %.2e\n', dataset_metrics.theta_cap_values(theta_idx));
    fprintf('Area Ratio: %.3f\n', dataset_metrics.area_ratio_values(area_idx));
    fprintf('SNR: %.1f\n', dataset_metrics.SNR_values(snr_idx));
    
    % Extract data for the selected point
    Y = dataset_metrics.Y{theta_idx, area_idx, snr_idx};
    Y_clean = dataset_metrics.Y_clean{theta_idx, area_idx, snr_idx};
    A0 = dataset_metrics.A0{theta_idx, area_idx, snr_idx};
    X0 = dataset_metrics.X0{theta_idx, area_idx, snr_idx};
    Aout = dataset_metrics.Aout{theta_idx, area_idx, snr_idx};
    Xout = dataset_metrics.Xout{theta_idx, area_idx, snr_idx};
    bout = dataset_metrics.bout{theta_idx, area_idx, snr_idx};
    extras = dataset_metrics.extras{theta_idx, area_idx, snr_idx};
    
    % Check if data exists for this point
    if isempty(Aout) || isempty(Xout)
        error('No reconstruction data available for this parameter combination');
    end
    
    % Add clean observation and demixing score to extras
    extras.Y_clean = Y_clean;
    [demix_score, overlap_matrix] = computeDemixingMetric(Xout);
    extras.demixing_score = demix_score;
    extras.demixing_matrix = overlap_matrix;
    
    % Create new figure for detailed visualization
    figure('Name', 'Detailed Results', 'Position', [100 100 1200 800]);
    
    % Call visualizeResults
    visualizeResults(Y, A0, Aout, X0, Xout, bout, extras);
    
    % Print performance metrics
    fprintf('\nPerformance Metrics:\n');
    fprintf('Kernel Similarity Score: %.3f\n', ...
        dataset_metrics.kernel_quality_final(theta_idx,area_idx,snr_idx));
    fprintf('Activation Score: %.3f\n', ...
        dataset_metrics.combined_activationScore(theta_idx,area_idx,snr_idx));
    fprintf('Demixing Score: %.4f\n', demix_score);
    
    % Print convergence information if available
    if isfield(extras, 'runtime')
        fprintf('\nConvergence Information:\n');
        fprintf('Runtime: %.2f seconds\n', extras.runtime);
    end
    
    if isfield(extras, 'phase1') && isfield(extras.phase1, 'residuals')
        fprintf('Final Residual: %.2e\n', extras.phase1.residuals(end));
    end
end

function [x,y,z] = ginput3d(n)
    % Get current axes
    ax = gca;
    
    % Get current view matrix
    [az,el] = view;
    
    % Get mouse click
    [x_2d,y_2d] = ginput(n);
    
    if isempty(x_2d)
        x = []; y = []; z = [];
        return;
    end
    
    % Get all scatter points
    h = findobj(ax, 'Type', 'scatter3');
    x_data = get(h, 'XData');
    y_data = get(h, 'YData');
    z_data = get(h, 'ZData');
    
    % Project all 3D points to 2D
    [x_proj,y_proj] = project_points(x_data,y_data,z_data,az,el);
    
    % Find nearest point
    dist = sqrt((x_proj-x_2d).^2 + (y_proj-y_2d).^2);
    [~,idx] = min(dist);
    
    % Return 3D coordinates of nearest point
    x = x_data(idx);
    y = y_data(idx);
    z = z_data(idx);
end

function [x_proj,y_proj] = project_points(x,y,z,az,el)
    % Project 3D points to 2D based on current view
    az = az * pi/180;
    el = el * pi/180;
    
    % Simple projection matrix
    x_proj = x.*cos(az) + y.*sin(az);
    y_proj = x.*sin(az).*sin(el) + y.*cos(az).*sin(el) + z.*cos(el);
end 