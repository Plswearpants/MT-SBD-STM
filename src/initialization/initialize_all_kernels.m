function A1_all = initialize_all_kernels(datasets)
% Initialize kernels for all datasets using consistent settings
%
% Input:
%   datasets: struct array containing synthetic datasets
%
% Output:
%   A1_all: cell array where each cell contains initialized kernels for one dataset
%
% Example:
%   A1_all = initialize_all_kernels(datasets);
%   A1_for_dataset_1 = A1_all{1};

fprintf('Initializing kernels for all datasets...\n');

% Get number of datasets
num_datasets = length(datasets);

% Initialize storage
A1_all = cell(num_datasets, 1);

% Set consistent initialization parameters
kerneltype = 'selected';
window_type = {'gaussian', 2};

% Initialize for each dataset
for i = 1:num_datasets
    fprintf('Dataset %d/%d: ', i, num_datasets);
    
    % Get dataset parameters
    Y = datasets(i).Y;
    kernel_sizes = datasets(i).params.kernel_sizes;
    num_kernels = size(kernel_sizes, 1);
    
    % Initialize kernels
    A1_all{i} = initialize_kernels(Y, num_kernels, kernel_sizes, kerneltype, window_type);
    
    fprintf('Done\n');
end

fprintf('All kernels initialized.\n');
end 