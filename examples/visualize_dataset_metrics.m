% Load the metrics
metrics = load_datasets_metrics();

%% Create the heatspace
metrics2heatspace(metrics);

%% view detailed dataset
visualize_heatspace_details(metrics);

%% load the metric with properGen results 
metrics= loadMetricDataset_new();