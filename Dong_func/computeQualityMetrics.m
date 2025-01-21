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

    % Set default indices if not provided
    if nargin < 6
        dataset_idx = [];
        param_idx = [];
    end

    % Create indices array for evaluateActivationReconstruction
    if ~isempty(dataset_idx) && ~isempty(param_idx)
        indices = [dataset_idx, param_idx];
    else
        indices = [];
    end

    % Compute activation metrics (with visualization off)
    activation_metric = evaluateActivationReconstruction(X0, Xout, kernel_size, false, indices);
    activation_similarity = activation_metric.similarity;
    
    % Compute kernel metrics (with visualization off)
    kernel_metric = evaluateKernelQuality(Aout, A0, false);
    kernel_similarity = kernel_metric;
end 