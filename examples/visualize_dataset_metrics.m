% Load the metrics
metrics = load_datasets_metrics();

%% load the metric with properGen results 
metrics= loadMetricDataset_new();

%% create heatspace
metrics2heat_properGen(metrics);

%% view detailed dataset
visualize_heatspace_details_properGen(metrics)

%% Create the heatspace
metrics2heatspace(metrics);

%% view detailed dataset
visualize_heatspace_details(metrics);

