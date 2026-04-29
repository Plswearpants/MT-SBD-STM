function visualize_heatspace_details_properGen(dataset_metrics, mode)
    if nargin < 2 || isempty(mode)
        mode = 1;
    end
    loader_mode = get_loader_axis_mode(dataset_metrics);

    % Show reference heatspace with the chosen mode
    metrics2heat_properGen(dataset_metrics, mode);
    
    has_repetitions = isfield(dataset_metrics, 'repetition_values') && ...
                      length(dataset_metrics.repetition_values) > 1;

    % Prompt user for parameter selection (axis-3 depends on mode)
    fprintf('\nAvailable parameter ranges:\n');
    fprintf('Theta Cap: [%.2e to %.2e]\n', min(dataset_metrics.theta_cap_values), max(dataset_metrics.theta_cap_values));
    fprintf('SNR: [%.1f to %.1f]\n', min(dataset_metrics.SNR_values), max(dataset_metrics.SNR_values));
    if mode == 2 && isfield(dataset_metrics, 'side_length_ratio_values') && ~isempty(dataset_metrics.side_length_ratio_values)
        fprintf('Side-length ratio: [%.4f to %.4f]\n', ...
            min(dataset_metrics.side_length_ratio_values), max(dataset_metrics.side_length_ratio_values));
    else
        fprintf('Nobs: [%.0f to %.0f]\n', min(dataset_metrics.Nobs_values), max(dataset_metrics.Nobs_values));
    end
    if has_repetitions
        fprintf('Repetitions: [%d to %d]\n', min(dataset_metrics.repetition_values), max(dataset_metrics.repetition_values));
    end
    
    theta_cap = input('Enter Theta Cap value: ');
    snr = input('Enter SNR value: ');

    if mode == 2 && isfield(dataset_metrics, 'side_length_ratio_values') && ~isempty(dataset_metrics.side_length_ratio_values)
        ratio_input = input('Enter Side-length ratio value: ');
        [~, side_idx] = min(abs(dataset_metrics.side_length_ratio_values - ratio_input));
        fprintf('  → Nearest ratio on grid: %.4f\n', dataset_metrics.side_length_ratio_values(side_idx));
    else
        ratio_input = [];
    end

    nobs_input = [];
    if isempty(ratio_input)
        nobs_input = input('Enter Nobs value: ');
    end
    
    if has_repetitions
        rep = input(sprintf('Enter Repetition number [%d-%d] (or press Enter to use first available): ', ...
            min(dataset_metrics.repetition_values), max(dataset_metrics.repetition_values)));
        if isempty(rep)
            rep = 1;
        end
    else
        rep = 1;
    end
    
    % Map inputs to indices in the primary [SNR × theta × axis3 × rep] arrays
    [~, theta_idx] = min(abs(dataset_metrics.theta_cap_values - theta_cap));
    [~, snr_idx] = min(abs(dataset_metrics.SNR_values - snr));

    axis_idx = [];
    if ~isempty(nobs_input)
        [~, axis_idx] = min(abs(dataset_metrics.Nobs_values - nobs_input));
    else
        if loader_mode == 2
            axis_idx = side_idx;
        else
            % Legacy mode-1 payload: side_length_ratio -> nearest populated N_obs
            axis_idx = resolve_nobs_from_side_ratio(dataset_metrics, snr_idx, theta_idx, side_idx, rep);
            if isempty(axis_idx)
                error('No reconstruction data found for the selected (SNR, density, side-length ratio) combination.');
            end
        end
    end
    
    if has_repetitions
        [~, rep_idx] = min(abs(dataset_metrics.repetition_values - rep));
        rep_actual = dataset_metrics.repetition_values(rep_idx);
    else
        rep_idx = 1;
        rep_actual = 1;
    end
    
    fprintf('\nUsing nearest available parameters:\n');
    fprintf('Theta Cap: %.2e\n', dataset_metrics.theta_cap_values(theta_idx));
    fprintf('SNR: %.1f\n', dataset_metrics.SNR_values(snr_idx));
    if loader_mode == 2 && isfield(dataset_metrics, 'Nobs_at_axis3')
        nobs_val = dataset_metrics.Nobs_at_axis3(snr_idx, theta_idx, axis_idx, rep_idx);
        fprintf('Nobs: %.0f\n', nobs_val);
    else
        fprintf('Nobs: %.0f\n', dataset_metrics.Nobs_values(axis_idx));
    end
    if ~isempty(ratio_input)
        fprintf('Side-length ratio (requested): %.4f\n', dataset_metrics.side_length_ratio_values(side_idx));
    end
    if has_repetitions
        fprintf('Repetition: %d\n', rep_actual);
    end
    
    % Extract data for the selected point (handle both 3D and 4D indexing)
    if has_repetitions && ndims(dataset_metrics.Y) == 4
        % 4D indexing: [SNR × theta × N_obs × rep]
        Y = dataset_metrics.Y{snr_idx, theta_idx, axis_idx, rep_idx};
        Y_clean = dataset_metrics.Y_clean{snr_idx, theta_idx, axis_idx, rep_idx};
        A0 = dataset_metrics.A0{snr_idx, theta_idx, axis_idx, rep_idx};
        A0_noiseless = dataset_metrics.A0_noiseless{snr_idx, theta_idx, axis_idx, rep_idx};
        X0 = dataset_metrics.X0{snr_idx, theta_idx, axis_idx, rep_idx};
        Aout = dataset_metrics.Aout{snr_idx, theta_idx, axis_idx, rep_idx};
        Xout = dataset_metrics.Xout{snr_idx, theta_idx, axis_idx, rep_idx};
        bout = dataset_metrics.bout{snr_idx, theta_idx, axis_idx, rep_idx};
        extras = dataset_metrics.extras{snr_idx, theta_idx, axis_idx, rep_idx};
    else
        % 3D indexing: [SNR × theta × N_obs] (backward compatibility)
        Y = dataset_metrics.Y{snr_idx, theta_idx, axis_idx};
        Y_clean = dataset_metrics.Y_clean{snr_idx, theta_idx, axis_idx};
        A0 = dataset_metrics.A0{snr_idx, theta_idx, axis_idx};
        A0_noiseless = dataset_metrics.A0_noiseless{snr_idx, theta_idx, axis_idx};
        X0 = dataset_metrics.X0{snr_idx, theta_idx, axis_idx};
        Aout = dataset_metrics.Aout{snr_idx, theta_idx, axis_idx};
        Xout = dataset_metrics.Xout{snr_idx, theta_idx, axis_idx};
        bout = dataset_metrics.bout{snr_idx, theta_idx, axis_idx};
        extras = dataset_metrics.extras{snr_idx, theta_idx, axis_idx};
    end
    
    % Check if data exists for this point
    if isempty(Aout) || isempty(Xout)
        error('No reconstruction data available for this parameter combination');
    end
    
    % Add clean observation and demixing score to extras
    extras.Y_clean = Y_clean;
    [demix_score, overlap_matrix] = computeDemixingMetric(Xout);
    extras.demixing_score = demix_score;
    extras.demixing_matrix = overlap_matrix;

    % Call visualizeResults
    visualizeResults(Y, A0_noiseless, Aout, X0, Xout, bout, extras);
    
    % Print performance metrics (handle both 3D and 4D indexing)
    fprintf('\nPerformance Metrics:\n');
    if has_repetitions && ndims(dataset_metrics.kernel_quality_final) == 4
        % Show metrics for this specific repetition
        fprintf('Kernel Quality Score(to noiseless GT): %.3f (Repetition %d)\n', ...
            dataset_metrics.kernel_quality_final(snr_idx, theta_idx, axis_idx, rep_idx), rep_actual);
        fprintf('Activation Similarity Score: %.3f (Repetition %d)\n', ...
            dataset_metrics.activation_similarity_final(snr_idx, theta_idx, axis_idx, rep_idx), rep_actual);
        fprintf('Combined Activation Score: %.3f (Repetition %d)\n', ...
            dataset_metrics.combined_activationScore(snr_idx, theta_idx, axis_idx, rep_idx), rep_actual);
        
        % Also show average over all repetitions
        rep_slice_kernel = dataset_metrics.kernel_quality_final(snr_idx, theta_idx, axis_idx, :);
        rep_slice_activation = dataset_metrics.activation_similarity_final(snr_idx, theta_idx, axis_idx, :);
        rep_slice_combined = dataset_metrics.combined_activationScore(snr_idx, theta_idx, axis_idx, :);
        avg_kernel = mean(rep_slice_kernel(:), 'omitnan');
        avg_activation = mean(rep_slice_activation(:), 'omitnan');
        avg_combined = mean(rep_slice_combined(:), 'omitnan');
        fprintf('\nAverage over all repetitions:\n');
        fprintf('  Kernel Quality Score: %.3f\n', avg_kernel);
        fprintf('  Activation Similarity Score: %.3f\n', avg_activation);
        fprintf('  Combined Activation Score: %.3f\n', avg_combined);
    else
        % 3D indexing (backward compatibility)
        fprintf('Kernel Quality Score(to noiseless GT): %.3f\n', ...
            dataset_metrics.kernel_quality_final(snr_idx, theta_idx, axis_idx));
        fprintf('Activation Similarity Score: %.3f\n', ...
            dataset_metrics.activation_similarity_final(snr_idx, theta_idx, axis_idx));
        fprintf('Combined Activation Score: %.3f\n', ...
            dataset_metrics.combined_activationScore(snr_idx, theta_idx, axis_idx));
    end
    fprintf('Demixing Score: %.4f\n', demix_score);
    
    % Print convergence information if available
    if isfield(extras, 'runtime')
        fprintf('\nConvergence Information:\n');
        fprintf('Runtime: %.2f seconds\n', extras.runtime);
    end
    
    if isfield(extras, 'phase1') && isfield(extras.phase1, 'residuals')
        fprintf('Final Residual: %.2e\n', extras.phase1.residuals(end));
    end
    
    % Print trajectory information if available
    if isfield(extras, 'phase1')
        fprintf('\nPhase I Trajectory Information:\n');
        if isfield(extras.phase1, 'kernel_quality_factors')
            kq_traj = extras.phase1.kernel_quality_factors;
            fprintf('Final Kernel Quality: %.3f\n', kq_traj(end));
        end
        if isfield(extras.phase1, 'activation_metrics')
            act_traj = extras.phase1.activation_metrics;
            fprintf('Final Activation Similarity: %.3f\n', act_traj(end));
        end
    end
end

function loader_mode = get_loader_axis_mode(dataset_metrics)
    if isfield(dataset_metrics, 'axis3_mode') && isscalar(dataset_metrics.axis3_mode)
        loader_mode = dataset_metrics.axis3_mode;
    else
        loader_mode = 1;
    end
end

function nobs_idx = resolve_nobs_from_side_ratio(dm, snr_idx, theta_idx, side_idx, rep_val)
%RESOLVE_NOBS_FROM_SIDE_RATIO Find N_obs index whose _by_side_length_ratio
%   slot matches the chosen side_idx and has reconstruction data.
    nobs_idx = [];
    has_rep = ndims(dm.Y) == 4;
    rep_idx = max(1, rep_val);

    % Build the set of candidate N_obs indices that have data
    for ni = 1:numel(dm.Nobs_values)
        if has_rep
            entry = dm.Aout{snr_idx, theta_idx, ni, rep_idx};
        else
            entry = dm.Aout{snr_idx, theta_idx, ni};
        end
        if isempty(entry)
            continue;
        end
        % Check if this N_obs maps to the requested side_idx
        if isfield(dm, 'kernel_quality_final_by_side_length_ratio')
            side_val = dm.kernel_quality_final_by_side_length_ratio(snr_idx, theta_idx, side_idx, rep_idx);
            if has_rep
                nobs_val = dm.kernel_quality_final(snr_idx, theta_idx, ni, rep_idx);
            else
                nobs_val = dm.kernel_quality_final(snr_idx, theta_idx, ni);
            end
            if ~isnan(side_val) && ~isnan(nobs_val) && abs(side_val - nobs_val) < 1e-10
                nobs_idx = ni;
                return;
            end
        end
    end

    % Fallback: return first non-empty N_obs for this (SNR, theta)
    for ni = 1:numel(dm.Nobs_values)
        if has_rep
            entry = dm.Aout{snr_idx, theta_idx, ni, rep_idx};
        else
            entry = dm.Aout{snr_idx, theta_idx, ni};
        end
        if ~isempty(entry)
            nobs_idx = ni;
            return;
        end
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