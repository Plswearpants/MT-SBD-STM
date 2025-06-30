%% ~~~~~~~~~~~~~ProperGen Results vis~~~~~~~~~~~~~~~~~~~~~
%% load the metric with properGen results 
metrics= loadMetricDataset_new();

%% create heatspace
metrics2heat_properGen(metrics);

%% view detailed dataset
visualize_heatspace_details_properGen(metrics)

%% view combined set 
% Load your 3 different experiment runs
metrics1 = loadMetricDataset_new();  % First run
metrics2 = loadMetricDataset_new();  % Second run  
metrics3 = loadMetricDataset_new();  % Third run

% Create cell array of metrics
dataset_metrics_array = {metrics1, metrics2, metrics3};

% Optional: provide names for each run
run_names = {'Experiment 1', 'Experiment 2', 'Experiment 3'};

% Plot all runs together using the new function
metrics2heat_multiple_runs(dataset_metrics_array, run_names);

%% ~~~~~~~~~~~~~Old Results vis~~~~~~~~~~~~~~~~~~~~~
%% Load the metrics
metrics = load_datasets_metrics();

%% Create the heatspace
metrics2heatspace(metrics);

%% view detailed dataset
visualize_heatspace_details(metrics);

