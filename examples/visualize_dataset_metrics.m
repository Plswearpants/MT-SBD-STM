%% ~~~~~~~~~~~~~ProperGen Results vis~~~~~~~~~~~~~~~~~~~~~
%% load the metric with properGen results 
metrics= loadMetricDataset_new();

%% create heatspace
metrics2heat_properGen(metrics);

%% fixed density vis
% For kernel similarity
metrics2heat_by_defect_density(metrics, 'kernel');

% For combined activation score
metrics2heat_by_defect_density(metrics, 'combined');
%%
metrics2heat_by_defect_density(metrics, 'multiplied');
%% view detailed dataset
visualize_heatspace_details_properGen(metrics)

%% number of defect occurrence vs kernel similarity
plot_defects_snr_kernel_similarity(metrics);
%%
plot_defects_vs_kernel_similarity_by_snr(metrics);
%% view combined set 
% Load your 3 different experiment runs
metrics1 = loadMetricDataset_new();  % First run
metrics2 = loadMetricDataset_new();  % Second run  
metrics3 = loadMetricDataset_new();  % Third run

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
metrics2heat_multiple_runs(dataset_metrics_array, run_names);

%% ~~~~~~~~~~~~~Old Results vis~~~~~~~~~~~~~~~~~~~~~
%% Load the metrics
metrics = load_datasets_metrics();

%% Create the heatspace
metrics2heatspace(metrics);

%% view detailed dataset
visualize_heatspace_details(metrics);

