%% 1. Initial Setup
% Load LDoS simulation data for kernel selection
load('example_data/LDoS_multi_1_defects_20260325_105344.mat');
LDoS_sim = LDoS_result;

% Display the 3D LDoS data for selection
fprintf('Displaying 3D LDoS simulation data...\n');
figure;
d3gridDisplay(LDoS_sim, 'dynamic');
title('LDoS Simulation Data - Use for Kernel Selection');

% Get user input for kernel slice selection
fprintf('\nPlease select slice indices (space-separated numbers): ');
input_str = input('', 's');
sliceidx = str2num(input_str);
num_kernels = length(sliceidx);

% Mask selected slices at image center (single-defect case).
% UI is only used to pick the masking radius per selected slice.
LDoS_sim = apply_center_defect_mask_selected_slices(LDoS_sim, sliceidx);

% Fixed parameters
fixed_params.p_scale = 3;                       % Resolution factor
fixed_params.N_single = N;                      % Input lattice size
fixed_params.num_kernels = num_kernels;         % Number of kernels

% Define parameter ranges for param_sets
SNR_values = [1, 9];                  % Different noise levels
defect_density_values = logspace(log10(1e-3), log10(1e-1), 15);     % activation densities with even multiplicative seperation (linear in log)
rep = 3;  % Number of repetitions (random activation patterns) per (N_obs, rho_d) combination

% Use side-length ratio as the generation handle and convert to N_obs upfront.
% side_length_ratio = kernel_side_pixels / observation_side_pixels = 2*M / N_obs.
side_length_ratio_values = linspace(0.05,0.35,13);

% Create handle-space parameter matrix
[S, D, R] = meshgrid(SNR_values, defect_density_values, side_length_ratio_values);
param_sets_handle = [S(:), D(:), R(:)];  % Each row: [SNR, defect_density, side_length_ratio]

% Convert handle-space params to the original N_obs-space param_sets
param_sets = zeros(size(param_sets_handle, 1), 3);  % [SNR, defect_density, N_obs]
for i = 1:size(param_sets_handle, 1)
    SNR_i = param_sets_handle(i, 1);
    rho_i = param_sets_handle(i, 2);
    side_length_ratio_i = param_sets_handle(i, 3);
    N_obs_i = convert_side_length_ratio_to_n_obs(side_length_ratio_i, SNR_i, ...
        LDoS_sim, sliceidx, fixed_params.N_single);
    param_sets(i, :) = [SNR_i, rho_i, N_obs_i];
end

fprintf('Parameter space setup complete:\n');
fprintf('- SNR values: %d points from %.2f to %.2f\n', length(SNR_values), min(SNR_values), max(SNR_values));
fprintf('- Defect density values: %d points from %.2e to %.2e\n', length(defect_density_values), min(defect_density_values), max(defect_density_values));
fprintf('- Side-length ratio values: %d points from %.3e to %.3e (2M/N_obs)\n', ...
    length(side_length_ratio_values), min(side_length_ratio_values), max(side_length_ratio_values));
fprintf('- Converted N_obs values: %d unique points from %d to %d\n', ...
    length(unique(param_sets(:,3))), min(param_sets(:,3)), max(param_sets(:,3)));
fprintf('- Repetitions (rep): %d per (N_obs, rho_d) combination\n', rep);
fprintf('Total parameter combinations: %d\n', size(param_sets, 1));
fprintf('Total datasets (with repetitions): %d\n', size(param_sets, 1) * rep);

% Display and confirm kernel selection
confirm_kernel_selection(LDoS_sim, sliceidx);

%% 2. Generate Base Activations (Step 1)
fprintf('\nStep 1: Generating base activations for different defect densities...\n');
base_activations = struct('X0', {}, 'defect_density', {}, 'N_obs', {}, 'repetition', {});

% Create figure for activation review
act_fig = figure('Name', 'Base Activation Pattern Review', ...
                'Position', [100 100 800 800]);

% Get unique combinations of defect_density and N_obs
unique_combinations = unique(param_sets(:,[2 3]), 'rows');  % [rho_d, N_obs]

% Generate rep activations for each unique combination
activation_counter = 1;
for i = 1:size(unique_combinations, 1)
    rho_d = unique_combinations(i,1);
    N_obs = unique_combinations(i,2);
    
    fprintf('Generating %d activations for rho_d = %.2e, N_obs = %d (%d/%d unique combinations)\n', ...
        rep, rho_d, N_obs, i, size(unique_combinations,1));
    
    % Generate rep activation patterns for this combination
    for r = 1:rep
        is_confirmed = false;
        while ~is_confirmed
            fprintf('  Generating repetition %d/%d...\n', r, rep);
            
            % Generate activation maps for this N_obs and density
            X0 = generate_activation_maps(N_obs, rho_d, fixed_params.p_scale, fixed_params.num_kernels);
            
            % Check for overlap between activations (generalized for any number of kernels)
            overlap_fraction = mean(sum(X0,3) > 1, 'all');
            has_overlap = overlap_fraction > 0;
            
            % If there's overlap, regenerate automatically
            if has_overlap
                fprintf('    Overlap detected (%.4f%%), regenerating...\n', 100*overlap_fraction);
                continue;  % Skip to next iteration to regenerate
            end
            
            % Display for review (only for first repetition, then auto-generate)
            if r == 1
                % Clear any previous confirmed flag
                if evalin('base', 'exist(''activation_confirmed'', ''var'')')
                    evalin('base', 'clear activation_confirmed');
                end
                
                clf(act_fig);
                
                % Create RGB image combining up to 3 activations for visualization
                rgb_activation = zeros([N_obs*fixed_params.p_scale N_obs*fixed_params.p_scale 3]);
                for ch = 1:min(3, fixed_params.num_kernels)
                    rgb_activation(:,:,ch) = X0(:,:,ch);
                end
                
                % Display combined activations
                imagesc(rgb_activation);
                if fixed_params.num_kernels == 2
                    title(sprintf('Combined Activations (\\rho_d=%.2e, N_{obs}=%d, Rep 1/%d)\nRed: Kernel 1, Green: Kernel 2\nYellow: Overlap', ...
                        rho_d, N_obs, rep), 'FontSize', 12);
                elseif fixed_params.num_kernels == 3
                    title(sprintf('Combined Activations (\\rho_d=%.2e, N_{obs}=%d, Rep 1/%d)\nRed: Kernel 1, Green: Kernel 2, Blue: Kernel 3', ...
                        rho_d, N_obs, rep), 'FontSize', 12);
                else
                    title(sprintf('Combined Activations (\\rho_d=%.2e, N_{obs}=%d, Rep 1/%d)...', ...
                        rho_d, N_obs, rep), 'FontSize', 12, 'Interpreter', 'tex');
                end
                axis square;
                
                % Add colorbar with custom labels
                colorbar('Ticks', [0, 0.5, 1], ...
                        'TickLabels', {'None', '', 'Active'});
                
                % Add UI controls
                regenerate_btn = uicontrol('Parent', act_fig, ...
                                         'Style', 'pushbutton', ...
                                         'String', 'Regenerate', ...
                                         'Position', [20 20 100 40], ...
                                         'Callback', @(src,event) regenerate_dataset());
                
                confirm_btn = uicontrol('Parent', act_fig, ...
                                       'Style', 'pushbutton', ...
                                       'String', 'Confirm', ...
                                       'Position', [140 20 100 40], ...
                                       'Callback', @(src,event) confirm_dataset());
                
                % Add text showing activation statistics for all kernels
                activation_rates = zeros(1, fixed_params.num_kernels);
                for k = 1:fixed_params.num_kernels
                    activation_rates(k) = 100*mean(X0(:,:,k), 'all');
                end
                stats_str = sprintf('Activation Rates:\n');
                for k = 1:fixed_params.num_kernels
                    stats_str = [stats_str, sprintf('Kernel %d: %.4f%%\n', k, activation_rates(k))];
                end
                stats_str = [stats_str, sprintf('Overlap: %.4f%%\n', 100*overlap_fraction)];
                stats_str = [stats_str, sprintf('Will generate %d total repetitions', rep)];
                
                uicontrol('Parent', act_fig, ...
                         'Style', 'text', ...
                         'String', stats_str, ...
                         'Position', [260 20 200 60+20*(fixed_params.num_kernels-2)], ...
                         'BackgroundColor', get(act_fig, 'Color'), ...
                         'HorizontalAlignment', 'left');
                
                drawnow;
                uiwait(act_fig);
                % Check if confirmed was set by the callback
                if evalin('base', 'exist(''activation_confirmed'', ''var'')')
                    is_confirmed = evalin('base', 'activation_confirmed');
                    evalin('base', 'clear activation_confirmed');  % Clear the variable
                else
                    % If regenerate was clicked, confirmed variable won't exist, so continue loop
                    is_confirmed = false;
                end
            else
                % For repetitions 2 through rep, auto-confirm (no UI)
                is_confirmed = true;
            end
        end
        
        % Store activation
        base_activations(activation_counter).X0 = X0;
        base_activations(activation_counter).defect_density = rho_d;
        base_activations(activation_counter).N_obs = N_obs;
        base_activations(activation_counter).repetition = r;
        activation_counter = activation_counter + 1;
    end
end
close(act_fig);

%% 3. Create Dataset Variations
fprintf('\nGenerating datasets for all parameter combinations...\n');
total_datasets = size(param_sets, 1) * rep;
final_datasets = cell(total_datasets, 1);
descriptions = cell(total_datasets, 1);

dataset_counter = 1;
for i = 1:size(param_sets, 1)
    % Get parameters for this iteration
    SNR = param_sets(i,1);
    rho_d = param_sets(i,2);
    N_obs = param_sets(i,3);
    
    % Calculate M for this SNR
    M = get_max_cutoff_over_slices(LDoS_sim, sliceidx, SNR, fixed_params.N_single);
    
    % Recover the designed side-length ratio from the handle-space grid
    designed_side_length_ratio = param_sets_handle(i, 3);
    % Calculate actual side-length ratio after rounding N_obs
    % (kernel-side / observation-side, pixel or lattice units are equivalent)
    side_length_ratio = 2 * M / N_obs;
    
    fprintf('Processing combination %d/%d: SNR=%.1f, rho_d=%.2e, N_obs=%d (designed_ratio=%.4f, actual_ratio=%.4f)\n', ...
        i, size(param_sets,1), SNR, rho_d, N_obs, designed_side_length_ratio, side_length_ratio);
    
    % Find all base activations matching this (rho_d, N_obs) combination
    base_indices = find([base_activations.defect_density] == rho_d & [base_activations.N_obs] == N_obs);
    
    if length(base_indices) ~= rep
        error('Expected %d activations for rho_d=%.2e, N_obs=%d, but found %d', ...
            rep, rho_d, N_obs, length(base_indices));
    end
    
    % Generate rep datasets, one for each activation pattern
    for r = 1:rep
        base_idx = base_indices(r);
        X0 = base_activations(base_idx).X0;
        repetition_idx = base_activations(base_idx).repetition;
        
        fprintf('  Generating dataset %d/%d (repetition %d/%d)...\n', ...
            dataset_counter, total_datasets, r, rep);
        
        % Generate kernels
        [A0, A0_noiseless] = generate_kernels(LDoS_sim(:,:,sliceidx), SNR, ...
            fixed_params.N_single, N_obs, fixed_params.p_scale);
        
        % Generate clean observation
        Y_clean = generate_clean_observation(A0_noiseless, X0);
        
        % Add noise
        [Y, A0] = add_noise_to_dataset(Y_clean, A0_noiseless, SNR);
        
        final_datasets{dataset_counter} = store_dataset(Y, Y_clean, A0, A0_noiseless, X0, randn, ...
            rho_d, SNR, N_obs, side_length_ratio, designed_side_length_ratio, repetition_idx);
        descriptions{dataset_counter} = sprintf('SNR=%.1f, ρ_d=%.2e, N_obs=%d, rep=%d (designed_ratio=%.4f, actual_ratio=%.4f)', ...
            SNR, rho_d, N_obs, repetition_idx, designed_side_length_ratio, side_length_ratio);
        dataset_counter = dataset_counter + 1;
    end
end

% Convert cell arrays to 1D structure arrays
datasets = [final_datasets{:}];
descriptions = reshape(descriptions, 1, []);

%% Save Results
% Add a note about the ordering in the saved file
ordering_info = struct('SNR_order', 'ascending', ...  % lowest to highest SNR
                      'SNR_values', SNR_values, ...
                      'rep', rep);  % Number of repetitions

% Save results with timestamp and ordering information
timestamp = datestr(now, 'yyyymmdd_HHMMSS');
save_dir = 'results/synthetic_datasets';
if ~exist(save_dir, 'dir')
    mkdir(save_dir);
end

save_filename = fullfile(save_dir, sprintf('synthetic_datasets_%s.mat', timestamp));
% Use -v7.3 (HDF5) so very large datasets (>2GB) can be saved.
save(save_filename, 'datasets', 'descriptions', 'param_sets', ...
    'fixed_params', 'sliceidx', 'LDoS_sim', 'ordering_info', 'rep', ...
    'side_length_ratio_values', '-v7.3');

fprintf('\nDatasets saved to: %s\n', save_filename);
fprintf('Saved %d datasets\n', numel(datasets));
fprintf('Parameter combinations: %d (SNR × defect_density × N_obs)\n', size(param_sets, 1));
fprintf('Repetitions per combination: %d\n', rep);
fprintf('Total datasets: %d\n', numel(datasets));

%% 6. Display Results Overview
fprintf('\nDataset Overview:\n');
fprintf('Total number of datasets: %d\n\n', numel(datasets));

% Get unique values for each parameter using arrayfun
unique_SNR = unique(arrayfun(@(x) x.params.SNR, datasets));
unique_rho_d = unique(arrayfun(@(x) x.params.defect_density, datasets));
unique_N_obs = unique(arrayfun(@(x) x.params.N_obs, datasets));
unique_repetitions = unique(arrayfun(@(x) x.params.repetition, datasets));

fprintf('Parameter ranges:\n');
fprintf('- SNR: %d values (%.1f to %.1f)\n', length(unique_SNR), min(unique_SNR), max(unique_SNR));
fprintf('- Defect density: %d values (%.2e to %.2e)\n', length(unique_rho_d), min(unique_rho_d), max(unique_rho_d));
fprintf('- N_obs: %d values (%d to %d)\n', length(unique_N_obs), min(unique_N_obs), max(unique_N_obs));
fprintf('- Repetitions: %d values (%d to %d)\n', length(unique_repetitions), min(unique_repetitions), max(unique_repetitions));

% Create separate figures for different N_obs values
for n = 1:length(unique_N_obs)
    % Get datasets for this N_obs
    idx_N_obs = find(arrayfun(@(x) x.params.N_obs == unique_N_obs(n), datasets));
    N_current = length(idx_N_obs);
    
    % Create 3D array for this N_obs
    obs_size = size(datasets(idx_N_obs(1)).Y);
    all_obs_N = zeros([obs_size N_current]);
    
    for i = 1:N_current
        all_obs_N(:,:,i) = datasets(idx_N_obs(i)).Y;
    end
    
    % Create new figure for this N_obs
    figure('Name', sprintf('Synthetic Observations N_obs=%d', unique_N_obs(n)));
    d3gridDisplay(all_obs_N, 'dynamic');
    title(sprintf('N_{obs} = %d (%d datasets)', unique_N_obs(n), N_current));
end

% Add a summary of side-length ratios (stored in params.area_ratio for compatibility)
figure('Name', 'Side-Length Ratio Distribution');
area_ratios = arrayfun(@(x) x.params.area_ratio, datasets);
histogram(area_ratios, 20);
title('Distribution of Side-Length Ratios');
xlabel('Side-Length Ratio (M/N_{obs})');
ylabel('Count');

% Print summary statistics
fprintf('\nSide-Length Ratio Statistics:\n');
fprintf('- Min: %.2e\n', min(area_ratios));
fprintf('- Max: %.2e\n', max(area_ratios));
fprintf('- Mean: %.2e\n', mean(area_ratios));
fprintf('- Median: %.2e\n', median(area_ratios));

%% Helper Functions
function X0 = generate_activation_maps(N_single, rho_d, p_scale, num_kernels)
    % Generate activation maps sequentially to ensure no overlap
    % Each kernel can only activate sites that are not already occupied
    
    % Initialize output
    X0 = zeros(N_single*p_scale, N_single*p_scale, num_kernels);
    
    % Calculate target number of defects for each kernel
    target_defects = round(N_single * N_single * rho_d);
    tolerance = 0.1;  % 10% tolerance
    hard_min_defects = 2;  % Hard lower bound per defect type
    computed_min_defects = max(round(target_defects * (1-tolerance)), 1);
    if computed_min_defects < hard_min_defects
        warning(['Computed min_defects=%d is below hard lower bound=%d. ' ...
            'Using hard lower bound for sampling.'], ...
            computed_min_defects, hard_min_defects);
    end
    min_defects = max(computed_min_defects, hard_min_defects);
    max_defects = max(round(target_defects * (1+tolerance)), min_defects);
    
    % Initialize occupancy map (tracks which sites are already occupied)
    occupancy_map = false(N_single*p_scale, N_single*p_scale);
    
    % Generate activations sequentially for each kernel
    for k = 1:num_kernels
        % Calculate how many defects we need for this kernel
        num_defects_needed = randi([min_defects, max_defects]);
        
        % Find all unoccupied sites
        [unoccupied_rows, unoccupied_cols] = find(~occupancy_map);
        
        if length(unoccupied_rows) < num_defects_needed
            warning('Not enough unoccupied sites for kernel %d. Need %d, have %d. Using all available.', ...
                k, num_defects_needed, length(unoccupied_rows));
            num_defects_needed = length(unoccupied_rows);
        end
        
        if num_defects_needed == 0
            warning('No unoccupied sites available for kernel %d', k);
            continue;
        end
        
        % Randomly select sites from unoccupied positions
        selected_indices = randperm(length(unoccupied_rows), num_defects_needed);
        selected_rows = unoccupied_rows(selected_indices);
        selected_cols = unoccupied_cols(selected_indices);
        
        % Place activations for this kernel
        for i = 1:num_defects_needed
            X0(selected_rows(i), selected_cols(i), k) = 1;
            occupancy_map(selected_rows(i), selected_cols(i)) = true;  % Mark as occupied
        end
        
        fprintf('Kernel %d: Placed %d activations on unoccupied sites\n', k, num_defects_needed);
    end
    
    % Verify no overlap (sanity check)
    total_activations = sum(X0, 3);
    max_overlap = max(total_activations(:));
    if max_overlap > 1
        error('Overlap detected! Maximum overlap is %d at some sites', max_overlap);
    else
        fprintf('Successfully generated non-overlapping activation maps for %d kernels. Max overlap: %d\n', num_kernels, max_overlap);
    end
end

function [A0, A0_noiseless] = generate_kernels(rho_single, SNR, N_single, N_obs, p_scale)
    % Generate kernels using SNR-based cutoff and proper scaling
    % 
    % Inputs:
    %   rho_single: Single defect QPI simulation (3D array for multiple slices)
    %   SNR: Signal-to-noise ratio
    %   N_single: Lattice size of input QPI pattern
    %   N_obs: Size of observation lattice
    %   p_scale: Resolution factor (pixels per lattice site)
    %
    % Outputs:
    %   A0: Cell array of noisy kernels
    %   A0_noiseless: Cell array of noiseless kernels
    
    num_kernels = size(rho_single, 3);
    A0 = cell(1, num_kernels);
    A0_noiseless = cell(1, num_kernels);
    
    % Process each kernel
    for k = 1:num_kernels
        % Find cutoff M for this slice
        [M, M_pixels, ~] = find_cutoff_noise_intersection(rho_single(:,:,k), SNR, N_single, false);
        
        % Get center and crop rho_single
        [ny, nx] = size(rho_single(:,:,k));
        center = ceil([ny, nx]/2);
        
        % Crop the pattern
        range_y = max(1, center(1)-M_pixels):min(ny, center(1)+M_pixels);
        range_x = max(1, center(2)-M_pixels):min(nx, center(2)+M_pixels);
        rho_single_cutoff = rho_single(range_y, range_x, k);
        
        % Resize cutoff pattern to match observation scale
        target_size = round([2*M*p_scale 2*M*p_scale]);
        A0_noiseless{k} = imresize(rho_single_cutoff, target_size);
        A0_noiseless{k} = proj2oblique(A0_noiseless{k});
        
        % Add noise based on SNR
        signal_variance = var(A0_noiseless{k}(:));
        eta = signal_variance / SNR;
        A0{k} = A0_noiseless{k} + sqrt(eta) * randn(size(A0_noiseless{k}));
        A0{k} = proj2oblique(A0{k});
    end
end

function [Y, A0] = add_noise_to_dataset(Y_clean, A0_noiseless, SNR)
    % Add noise to clean observation and kernels
    % Initialize noisy kernels
    A0 = cell(size(A0_noiseless));

    % Auto-estimate signal scale from center cuts (0 deg and 45 deg) per kernel.
    % For each cut, use first peak prominence as peak-to-valley amplitude.
    p2v_per_kernel = compute_p2v_centercuts_per_kernel(A0_noiseless);
    mean_p2v = mean(p2v_per_kernel, 'omitnan');
    if ~isfinite(mean_p2v) || mean_p2v <= 0
        error('Failed to compute positive mean peak-to-valley amplitude from kernels.');
    end

    % Keep SNR convention aligned with real-data pipeline:
    % SNR = mean_p2v / sqrt(eta)  =>  eta = (mean_p2v / SNR)^2
    eta = (mean_p2v / SNR)^2;
    
    % Add noise to kernels
    for k = 1:length(A0_noiseless)
        A0{k} = A0_noiseless{k} + sqrt(eta) * randn(size(A0_noiseless{k}));
        A0{k} = proj2oblique(A0{k});
    end
    
    % Add noise to observation using the same eta from p2v/SNR rule.
    Y = Y_clean + sqrt(eta) * randn(size(Y_clean));
end

function p2v_per_kernel = compute_p2v_centercuts_per_kernel(A0_noiseless)
%COMPUTE_P2V_CENTERCUTS_PER_KERNEL Estimate p2v amplitude per kernel.
%   For each kernel:
%   - take horizontal (0 deg) and diagonal (45 deg) center cuts
%   - for each cut, pick the highest-prominence peak
%   - use that prominence as peak-to-valley amplitude
%   - take max of the two cuts to get one p2v per kernel

    num_kernels = numel(A0_noiseless);
    p2v_per_kernel = nan(1, num_kernels);

    for k = 1:num_kernels
        Ak = A0_noiseless{k};
        cut0 = center_cut_profile(Ak, 0);   % 0 deg center cut
        cut45 = center_cut_profile(Ak, 45); % 45 deg center cut

        p2v0 = first_peak_prominence(cut0);
        p2v45 = first_peak_prominence(cut45);

        p2v_per_kernel(k) = max([p2v0, p2v45], [], 'omitnan');
    end
end

function prof = center_cut_profile(A, angle_deg)
%CENTER_CUT_PROFILE Sample a center-passing line profile at given angle.
    [h, w] = size(A);
    cx = (w + 1) / 2;
    cy = (h + 1) / 2;
    dx = cosd(angle_deg);
    dy = sind(angle_deg);

    % Symmetric half-length that stays inside bounds for +/- direction.
    lims = [];
    if abs(dx) > eps
        lims(end+1) = (w - cx) / abs(dx); %#ok<AGROW>
        lims(end+1) = (cx - 1) / abs(dx); %#ok<AGROW>
    end
    if abs(dy) > eps
        lims(end+1) = (h - cy) / abs(dy); %#ok<AGROW>
        lims(end+1) = (cy - 1) / abs(dy); %#ok<AGROW>
    end
    half_len = max(1, min(lims));

    x1 = cx - half_len * dx;
    y1 = cy - half_len * dy;
    x2 = cx + half_len * dx;
    y2 = cy + half_len * dy;

    ns = max(2, round(2 * half_len) + 1);
    prof = improfile(A, [x1 x2], [y1 y2], ns, 'bilinear');
    prof = double(prof(:)).';
    prof = prof(~isnan(prof));
end

function p2v = first_peak_prominence(sig)
%FIRST_PEAK_PROMINENCE Return maximum peak prominence.
    sig = sig(:).';
    x = 1:numel(sig);

    [pks, ~, ~, proms] = findpeaks(sig, x, ...
        'Annotate', 'extents', ...
        'WidthReference', 'halfheight');

    if isempty(pks)
        % Fallback for monotonic/noisy edge cases: peak minus global valley.
        p2v = max(sig) - min(sig);
        return;
    end

    p2v = max(proms);
end

function Y_clean = generate_clean_observation(A0_noiseless, X0)
    % Generate clean observation by convolving kernels with activation maps
    Y_clean = zeros(size(X0, 1:2));
    
    for k = 1:size(X0, 3)
        Y_clean = Y_clean + convfft2(A0_noiseless{k}, X0(:,:,k));
    end
end

function confirm_kernel_selection(LDoS_sim, sliceidx)
    figure('Name', 'Selected Kernel Slices');
    for k = 1:length(sliceidx)
        subplot(1, length(sliceidx), k);
        imagesc(LDoS_sim(:,:,sliceidx(k)));
        title(sprintf('Kernel %d (Slice %d)', k, sliceidx(k)));
        colorbar;
        axis square;
    end
    
    confirmation = input('\nAre you satisfied with these kernel selections? (y/n): ', 's');
    while ~strcmpi(confirmation, 'y')
        % Allow reselection if not satisfied
        fprintf('\nPlease reselect slices:\n');
        for k = 1:length(sliceidx)
            sliceidx(k) = input(sprintf('Enter slice index for kernel %d: ', k));
            while sliceidx(k) < 1 || sliceidx(k) > size(LDoS_sim,3)
                fprintf('Invalid slice index. Please enter a number between 1 and %d\n', size(LDoS_sim,3));
                sliceidx(k) = input(sprintf('Enter slice index for kernel %d: ', k));
            end
        end
        
        % Update display
        for k = 1:length(sliceidx)
            subplot(1, length(sliceidx), k);
            imagesc(LDoS_sim(:,:,sliceidx(k)));
            title(sprintf('Kernel %d (Slice %d)', k, sliceidx(k)));
            colorbar;
            axis square;
        end
        
        confirmation = input('\nAre you satisfied with these kernel selections? (y/n): ', 's');
    end
end

function N_obs = convert_side_length_ratio_to_n_obs(side_length_ratio, SNR, LDoS_sim, sliceidx, N_single)
    % Convert (side_length_ratio, SNR, selected kernel slices) to N_obs
    % using the largest cutoff over all selected slices.
    M_max = get_max_cutoff_over_slices(LDoS_sim, sliceidx, SNR, N_single);
    N_obs = max(3, round((2 * M_max) / side_length_ratio));
end

function M_max = get_max_cutoff_over_slices(LDoS_sim, sliceidx, SNR, N_single)
    M_values = zeros(1, length(sliceidx));
    for k = 1:length(sliceidx)
        [M_values(k), ~] = find_cutoff_noise_intersection(LDoS_sim(:,:,sliceidx(k)), SNR, N_single, false);
    end
    M_max = max(M_values);
end

function dataset = store_dataset(Y, Y_clean, A0, A0_noiseless, X0, b0, rho_d, SNR, N_obs, area_ratio, designed_side_length_ratio, repetition)
    dataset.Y = Y;
    dataset.Y_clean = Y_clean;
    dataset.A0 = A0;
    dataset.A0_noiseless = A0_noiseless;
    dataset.X0 = X0;
    
    % Calculate kernel sizes as 2×2 matrix where each row is [height, width] for one kernel
    kernel_sizes = zeros(2, 2);  % Initialize [2×2] matrix
    for k = 1:length(A0)
        kernel_sizes(k,:) = size(A0{k});  % Store [height, width] for each kernel
    end
    
    % Initialize A1 using most isolated points method
    dataset.A1 = initialize_kernels_from_isolated_points(Y, X0, A0);
    
    dataset.params = struct('defect_density', rho_d, ...
                          'kernel_size', kernel_sizes, ...
                          'SNR', SNR, ...
                          'N_obs', N_obs, ...
                          'area_ratio', area_ratio, ...
                          'designed_side_length_ratio', designed_side_length_ratio, ...
                          'repetition', repetition);
    dataset.b0 = b0;
end

function A1 = initialize_kernels_from_isolated_points(Y, X0, A0)
    % Initialize kernels by finding most isolated activation points
    % Inputs:
    %   Y: Observation data
    %   X0: Ground truth activation maps (binary)
    %   A0: Ground truth kernels (for size reference)
    % Output:
    %   A1: Cell array of initialized kernels with same sizes as A0
    
    num_kernels = size(X0, 3);
    A1 = cell(1, num_kernels);
    
    % Get kernel sizes from A0
    kernel_sizes = zeros(num_kernels, 2);
    for k = 1:num_kernels
        kernel_sizes(k,:) = size(A0{k});
    end
    
    % Find most isolated points for each kernel
    most_isolated_points = cell(1, num_kernels);
    
    for k = 1:num_kernels
        % Get positions of defects (all activations are 1 in synthetic data)
        [rows, cols] = find(X0(:,:,k));
        defect_positions = [rows, cols];
        
        if isempty(defect_positions)
            warning('No defects found for kernel %d', k);
            A1{k} = zeros(kernel_sizes(k,:));  % Initialize with zeros of correct size
            continue;
        end
        
        % Create summed activation map of all other kernels
        X_others = zeros(size(X0(:,:,1)));
        for l = 1:num_kernels
            if l ~= k
                X_others = X_others + X0(:,:,l);
            end
        end
        
        % Get positions of defects in other kernels
        [other_rows, other_cols] = find(X_others);
        other_positions = [other_rows, other_cols];
        
        % Calculate isolation scores for all points
        S_k = zeros(size(defect_positions, 1), 1);
        for i = 1:size(defect_positions, 1)
            diffs = other_positions - defect_positions(i,:);
            distances = sum(diffs.^2, 2);
            S_k(i) = min(distances);
        end
        
        % Get exact kernel size from A0
        [kernel_h, kernel_w] = size(A0{k});
        half_h = floor(kernel_h/2);
        half_w = floor(kernel_w/2);
        
        % Filter points that are away from boundaries
        valid_points = true(size(defect_positions, 1), 1);
        for i = 1:size(defect_positions, 1)
            y = defect_positions(i,1);
            x = defect_positions(i,2);
            
            if y <= half_h || y >= size(X0,1) - half_h || ...
               x <= half_w || x >= size(X0,2) - half_w
                valid_points(i) = false;
            end
        end
        
        % If we have valid points away from boundaries, use the one with highest score
        if any(valid_points)
            valid_defects = defect_positions(valid_points,:);
            valid_scores = S_k(valid_points);
            [~, max_idx] = max(valid_scores);
            most_isolated_points{k} = valid_defects(max_idx,:);
            
            % Extract kernel directly since we know it's away from boundaries
            y = most_isolated_points{k}(1);
            x = most_isolated_points{k}(2);
            A1{k} = Y(y-half_h:y+half_h, x-half_w:x+half_w);
            
            % Ensure exact size match with A0
            if ~isequal(size(A1{k}), size(A0{k}))
                A1{k} = imresize(A1{k}, size(A0{k}));
            end
        else
            % Fallback: use the point with highest isolation score and pad
            [~, max_idx] = max(S_k);
            most_isolated_points{k} = defect_positions(max_idx,:);
            
            % Extract kernel with padding
            y = most_isolated_points{k}(1);
            x = most_isolated_points{k}(2);
            
            % Calculate valid ranges that stay within image boundaries
            y_start = max(1, y - half_h);
            y_end = min(size(Y,1), y + half_h);
            x_start = max(1, x - half_w);
            x_end = min(size(Y,2), x + half_w);
            
            % Extract the valid portion
            kernel_patch = Y(y_start:y_end, x_start:x_end);
            
            % Create zero-padded kernel of required size
            A1{k} = zeros(kernel_sizes(k,:));
            
            % Calculate offsets for centering the patch
            y_offset = half_h - (y - y_start);
            x_offset = half_w - (x - x_start);
            
            % Place the patch in the center of the zero-padded kernel
            A1{k}(y_offset + (1:size(kernel_patch,1)), x_offset + (1:size(kernel_patch,2))) = kernel_patch;
            
            % Ensure exact size match with A0
            if ~isequal(size(A1{k}), size(A0{k}))
                A1{k} = imresize(A1{k}, size(A0{k}));
            end
        end
        
        % Normalize the kernel
        A1{k} = proj2oblique(A1{k});
        
        % Final size check
        assert(isequal(size(A1{k}), size(A0{k})), ...
            'Kernel size mismatch: A1{%d} size %dx%d does not match A0{%d} size %dx%d', ...
            k, size(A1{k},1), size(A1{k},2), k, size(A0{k},1), size(A0{k},2));
    end
end

function regenerate_dataset()
    uiresume;
end

function confirm_dataset()
    assignin('base', 'activation_confirmed', true);
    uiresume;
end 

function Y_out = apply_center_defect_mask_selected_slices(Y_in, sliceidx)
%APPLY_CENTER_DEFECT_MASK_SELECTED_SLICES Apply center Gaussian mask to chosen slices.
%   For each selected slice, user chooses only the radius using a fixed-center circle.

    Y_out = Y_in;
    [h, w, ~] = size(Y_in);
    center = [(w + 1) / 2, (h + 1) / 2]; % [x, y]

    fprintf('\nCenter-defect masking for selected slices...\n');
    for i = 1:numel(sliceidx)
        s = sliceidx(i);
        fprintf('  Slice %d (%d/%d): choose mask radius and confirm.\n', s, i, numel(sliceidx));

        sigma = select_center_radius_ui(Y_in(:,:,s), center);
        center_xy = center;
        mask2d = build_center_gaussian_mask([h, w], center_xy, sigma);

        % Use gaussianMaskDefects with provided mask (no center selection UI).
        slice_3d = reshape(Y_in(:,:,s), [h, w, 1]);
        [masked_slice_3d, ~] = gaussianMaskDefects(slice_3d, 1, 1, mask2d);
        Y_out(:,:,s) = masked_slice_3d(:,:,1);
    end
end

function sigma = select_center_radius_ui(slice_data, center)
%SELECT_CENTER_RADIUS_UI Radius-only UI with fixed center.

    f = figure('Name', 'Center Defect Mask Radius', 'Position', [100, 100, 800, 650]);
    imagesc(slice_data);
    colormap('gray'); colorbar; axis image;
    caxis([min(slice_data(:)), max(slice_data(:))]);
    title('Adjust radius at fixed center, then click Confirm');
    hold on;
    plot(center(1), center(2), 'r+', 'MarkerSize', 12, 'LineWidth', 1.5);

    default_radius = max(2, 0.08 * min(size(slice_data,1), size(slice_data,2)));
    h_circle = drawcircle('Center', center, 'Radius', default_radius, 'Color', 'r', 'FaceAlpha', 0.1);
    addlistener(h_circle, 'ROIMoved', @(src,~) keep_center_fixed(src, center));

    uicontrol('Style', 'text', 'String', 'Radius:', ...
              'Position', [10, 10, 60, 20], 'BackgroundColor', get(f,'Color'));
    h_edit = uicontrol('Style', 'edit', 'String', num2str(default_radius, '%.2f'), ...
                       'Position', [70, 10, 70, 24], ...
                       'Callback', @(src,~) update_radius_from_edit(src, h_circle));
    uicontrol('Style', 'pushbutton', 'String', 'Confirm', ...
              'Position', [150, 10, 90, 24], ...
              'Callback', @(~,~) uiresume(f));

    addlistener(h_circle, 'ROIMoved', @(src,~) sync_radius_edit(src, h_edit));
    uiwait(f);
    sigma = h_circle.Radius;
    if isvalid(f)
        close(f);
    end
end

function keep_center_fixed(h_circle, center)
    if isvalid(h_circle)
        h_circle.Center = center;
    end
end

function update_radius_from_edit(h_edit, h_circle)
    r = str2double(h_edit.String);
    if isfinite(r) && r > 0 && isvalid(h_circle)
        h_circle.Radius = r;
    end
end

function sync_radius_edit(h_circle, h_edit)
    if isvalid(h_circle) && isvalid(h_edit)
        h_edit.String = num2str(h_circle.Radius, '%.2f');
    end
end

function mask2d = build_center_gaussian_mask(sz_hw, center_xy, sigma)
%BUILD_CENTER_GAUSSIAN_MASK Build mask matching gaussianMaskDefects profile.

    h = sz_hw(1);
    w = sz_hw(2);
    [X, Y] = meshgrid(1:w, 1:h);
    distance_squared = (X - center_xy(1)).^2 + (Y - center_xy(2)).^2;
    distance = sqrt(distance_squared);

    step_loc = 2 * sigma;
    step_shapeness = 10;
    smooth_step = 0.5 + 0.5 * tanh(-step_shapeness * (distance - step_loc));
    gaussian = 0.99 * exp(-distance_squared / (2 * sigma^2)) .* smooth_step;
    mask2d = 1 - gaussian;
end