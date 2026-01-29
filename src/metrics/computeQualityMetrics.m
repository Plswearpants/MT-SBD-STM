function [activation_similarity, kernel_similarity] = computeQualityMetrics(X0, Xout, A0, Aout, kernel_size, dataset_idx, param_idx)
    % Bundles both activation and kernel quality metrics into one function
    % without visualization
    %
    % Inputs:
    %   X0: ground truth activation maps [height x width x num_kernels]
    %   Xout: reconstructed activation maps [height x width x num_kernels]
    %   A0: ground truth kernels (cell array)
    %   Aout: reconstructed kernels (cell array)
    %   kernel_size: size of each kernel [num_kernels x 2]
    %   dataset_idx: (optional) dataset index for labeling
    %   param_idx: (optional) parameter set index for labeling
    %
    % Outputs:
    %   activation_similarity: array of similarity scores for activations
    %   kernel_similarity: array of similarity scores for kernels

    % Input validation
    validateattributes(X0, {'numeric'}, {'3d'}, mfilename, 'X0');
    validateattributes(Xout, {'numeric'}, {'3d'}, mfilename, 'Xout');
    validateattributes(A0, {'cell'}, {'vector'}, mfilename, 'A0');
    validateattributes(Aout, {'cell'}, {'vector'}, mfilename, 'Aout');
    validateattributes(kernel_size, {'numeric'}, {'2d'}, mfilename, 'kernel_size');

    % Size consistency checks
    if size(X0, 3) ~= size(Xout, 3)
        error('Number of kernels mismatch between X0 and Xout');
    end
    if length(A0) ~= length(Aout)
        error('Number of kernels mismatch between A0 and Aout');
    end
    if size(kernel_size, 1) ~= length(A0) || size(kernel_size, 2) ~= 2
        error('kernel_size must be [num_kernels x 2] matrix');
    end

    % Set default indices if not provided
    if nargin < 6
        dataset_idx = [];
        param_idx = [];
    end

    % Validate optional indices when provided
    if ~isempty(dataset_idx)
        validateattributes(dataset_idx, {'numeric'}, {'scalar', 'positive'}, mfilename, 'dataset_idx');
    end
    if ~isempty(param_idx)
        validateattributes(param_idx, {'numeric'}, {'scalar', 'positive'}, mfilename, 'param_idx');
    end

    % Create indices array for evaluateActivationReconstruction
    if ~isempty(dataset_idx) && ~isempty(param_idx)
        indices = [dataset_idx, param_idx];
    else
        indices = [];
    end

    % Compute activation metrics (with visualization off)
    try
        activation_metric = evaluateActivationReconstruction(X0, Xout, kernel_size, false, indices);
        activation_similarity = activation_metric.similarity;
    catch ME
        error('Failed to compute activation metrics: %s', ME.message);
    end
    
    % Compute kernel metrics (with visualization off)
    try
        kernel_metric = evaluateKernelQuality(Aout, A0, false);
        kernel_similarity = kernel_metric;
    catch ME
        error('Failed to compute kernel metrics: %s', ME.message);
    end

    % Validate outputs
    if any(isnan(activation_similarity)) || any(isnan(kernel_similarity))
        warning('Quality metrics contain NaN values');
    end
    if any(activation_similarity < 0) || any(activation_similarity > 1) || ...
       any(kernel_similarity < 0) || any(kernel_similarity > 1)
        warning('Quality metrics outside expected range [0,1]');
    end
end 