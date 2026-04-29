%% ~~~~~~~~~~~~~ProperGen Results vis~~~~~~~~~~~~~~~~~~~~~
% Global axis-3 mode: 1 = N_obs, 2 = side_length_ratio.
% Set once here; all downstream visualizations follow this choice.
axis3_mode = 2; 
%% load the metric with properGen results 
metrics = loadMetricDataset_new(axis3_mode);
%% Build the observation fidelity 
metrics =build_observation_fidelity_metrics(metrics, axis3_mode);

%% Build the normalized kernel similarity 
metrics = build_normalized_kernel_similarity(metrics);

%% Generalized heatspace plot
metrics2heat_general(metrics, 2);

%% Line profile draw on metric heatspace
metrics2heat_general(metrics, 2, ...
    'plot_mode', 'line_profile', ...
    'metric_type', 'combined', ...
    'snr_value', 5, ...
    'interp_factor', 1, ...
    'manual_colormap', slanCM('viridis'));
%% create heatspace
metrics2heat_properGen(metrics, axis3_mode);

%% fixed density vis
% For kernel similarity
metrics2heat_by_defect_density(metrics, 'kernel', axis3_mode);
%% per-SNR 2D heatmaps with linear interpolation (same colormap)
metrics2heat_by_snr_interpolated( ...
    metrics, 'kernel', axis3_mode, 5);
%% per-SNR baseline heatmaps for kernel similarity normalization
metric_colormap = slanCM('viridis');
contour_levels =[0.95,0.85];

metrics2heat_by_snr_interpolated( ...
    metrics, 'kernel', axis3_mode, 5, metric_colormap, contour_levels);
%%
metrics2heat_by_snr_interpolated( ...
    metrics, 'combined', axis3_mode, 5, metric_colormap, contour_levels);
%% interactive explorer on fixed SNR heatmap (click to inspect nearest dataset)
explore_snr_heatmap_click(metrics, 5, 'combined', axis3_mode, 1);

%% For combined activation score
metrics2heat_by_defect_density(metrics, 'combined', axis3_mode);
%%
metrics2heat_by_defect_density(metrics, 'multiplied', axis3_mode);
%% view detailed dataset
visualize_heatspace_details_properGen(metrics, axis3_mode)

%% number of defect occurrence vs kernel similarity
plot_defects_snr_kernel_similarity(metrics);
%% number of defect occurrence vs kernel similarity by SNR
plot_defects_vs_kernel_similarity_by_snr(metrics);
%% given SNR: occurrence vs side-length ratio (2D scatter)
plot_occurrence_vs_length_ratio_by_snr(metrics);
%% Build NOR and LOO vs density metric (filter by SNR and side-length ratio)
metrics = build_nor_loo_metrics(metrics);
%% Plot NOR and LOO vs density on given values
selected_snr = [5];
selected_side_length_ratio = [0.175];
[nor_density_axis, nor_curve, loo_curve, nor_loo_curve] = plot_nor_loo_vs_density(metrics, selected_snr, selected_side_length_ratio);
%% view combined set 
% Load your 3 different experiment runs
metrics1 = loadMetricDataset_new(1);  % First run (combine_metrics_for_plotting expects N_obs axis)
metrics2 = loadMetricDataset_new(1);  % Second run  
metrics3 = loadMetricDataset_new(1);  % Third run

% Create cell array of metrics
dataset_metrics_array = {metrics1, metrics2, metrics3};

% Optional: provide names for each run
run_names = {'Experiment 1', 'Experiment 2', 'Experiment 3'};
%% number of defect occurrence vs kernel similarity(3D) -> combined version
combined_metrics = combine_metrics_for_plotting(dataset_metrics_array);
% Plot
plot_defects_snr_kernel_similarity(combined_metrics);
%% number of defect occurrence vs kernel similarity(2D)
combined_metrics = combine_metrics_for_plotting(dataset_metrics_array);

plot_defects_vs_kernel_similarity_by_snr(combined_metrics);
%% Plot all runs together using the new function
metrics2heat_multiple_runs(dataset_metrics_array, run_names, axis3_mode);

%% ~~~~~~~~~~~~~Old Results vis~~~~~~~~~~~~~~~~~~~~~
%% Load the metrics
metrics = load_datasets_metrics();

%% Create the heatspace
metrics2heatspace(metrics);

%% view detailed dataset
visualize_heatspace_details(metrics);

