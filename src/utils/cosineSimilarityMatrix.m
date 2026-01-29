function similarity_matrix = cosineSimilarityMatrix(matrices, varargin)
    % COSINESIMILARITYMATRIX Calculate cosine similarity matrix across n matrices
    %
    % This function computes the cosine similarity between all pairs of matrices
    % using evaluateActivationReconstruction.m as the underlying calculation method.
    % The result is a symmetric matrix where each entry (i,j) represents the 
    % cosine similarity between matrices i and j.
    %
    % Inputs:
    %   matrices: cell array of matrices to compare, or 3D array where the
    %             third dimension represents different matrices
    %   varargin: optional parameters
    %       'kernel_size': size of each kernel [num_matrices x 2] (required for 3D arrays)
    %       'visualize': boolean for visualization (default: false)
    %       'indices': [dataset_num, param_num] for labeling plots (optional)
    %
    % Outputs:
    %   similarity_matrix: n x n symmetric matrix of cosine similarities
    %                      where n is the number of input matrices
    %
    % Example:
    %   % Using cell array of matrices
    %   matrices = {rand(10,10), rand(10,10), rand(10,10)};
    %   sim_matrix = cosineSimilarityMatrix(matrices);
    %
    %   % Using 3D array
    %   matrices_3d = rand(10, 10, 3);
    %   kernel_size = [10, 10; 10, 10; 10, 10];  % Required for 3D arrays
    %   sim_matrix = cosineSimilarityMatrix(matrices_3d, 'kernel_size', kernel_size);
    %
    % Notes:
    %   - Uses evaluateActivationReconstruction.m for similarity computation
    %   - Includes alignment and density-adaptive Gaussian filtering
    %   - The diagonal elements will always be 1 (self-similarity)
    %   - The matrix is symmetric (sim_matrix(i,j) = sim_matrix(j,i))
    
    % Parse optional parameters
    p = inputParser;
    addParameter(p, 'kernel_size', [], @isnumeric);
    addParameter(p, 'visualize', false, @islogical);
    addParameter(p, 'indices', [], @isnumeric);
    parse(p, varargin{:});
    
    kernel_size = p.Results.kernel_size;
    visualize = p.Results.visualize;
    indices = p.Results.indices;
    
    % Input validation and preprocessing
    if iscell(matrices)
        % Handle cell array input
        num_matrices = length(matrices);
        
        % Validate all matrices have the same size
        if num_matrices < 2
            error('At least 2 matrices are required for similarity computation');
        end
        
        % Get size of first matrix for validation
        first_size = size(matrices{1});
        
        % Validate all matrices
        for i = 1:num_matrices
            if ~isnumeric(matrices{i})
                error('All matrices must be numeric');
            end
            
            if ~isequal(size(matrices{i}), first_size)
                error('All matrices must have the same size');
            end
        end
        
        % Convert cell array to 3D array for evaluateActivationReconstruction
        [rows, cols] = size(matrices{1});
        matrices_3d = zeros(rows, cols, num_matrices);
        for i = 1:num_matrices
            matrices_3d(:,:,i) = matrices{i};
        end
        
    elseif isnumeric(matrices) && ndims(matrices) == 3
        % Handle 3D array input
        [rows, cols, num_matrices] = size(matrices);
        matrices_3d = matrices;
        
        if num_matrices < 2
            error('At least 2 matrices are required for similarity computation');
        end
        
        % Validate kernel_size is provided for 3D arrays
        if isempty(kernel_size)
            error('kernel_size must be provided for 3D array input');
        end
        
        if size(kernel_size, 1) ~= num_matrices || size(kernel_size, 2) ~= 2
            error('kernel_size must be [num_matrices x 2] matrix');
        end
        
    else
        error('Input must be either a cell array of matrices or a 3D numeric array');
    end
    
    % If kernel_size not provided for cell array, use matrix size
    if isempty(kernel_size)
        kernel_size = repmat([rows, cols], num_matrices, 1);
    end
    
    % Initialize similarity matrix
    similarity_matrix = zeros(num_matrices, num_matrices);
    
    % Compute similarity for all pairs using evaluateActivationReconstruction
    for i = 1:num_matrices
        for j = i:num_matrices  % Only compute upper triangle + diagonal
            if i == j
                % Self-similarity is always 1
                similarity_matrix(i, j) = 1;
            else
                % Use evaluateActivationReconstruction for pair-wise similarity
                X0 = matrices_3d(:,:,i);
                Xout = matrices_3d(:,:,j);
                
                % Create 3D arrays for evaluateActivationReconstruction
                X0_3d = reshape(X0, [size(X0), 1]);
                Xout_3d = reshape(Xout, [size(Xout), 1]);
                
                % Use single kernel size for this pair
                pair_kernel_size = kernel_size(i, :);
                
                % Call evaluateActivationReconstruction
                if ~isempty(indices)
                    [metrics, ~] = evaluateActivationReconstruction(X0_3d, Xout_3d, pair_kernel_size, visualize, indices);
                else
                    [metrics, ~] = evaluateActivationReconstruction(X0_3d, Xout_3d, pair_kernel_size, visualize);
                end
                
                % Extract similarity score
                similarity_matrix(i, j) = metrics.similarity;
                
                % Fill in symmetric part
                similarity_matrix(j, i) = similarity_matrix(i, j);
            end
        end
    end
    
    % Validate output
    if any(isnan(similarity_matrix(:)))
        warning('Similarity matrix contains NaN values');
    end
    
    if any(similarity_matrix(:) < -1) || any(similarity_matrix(:) > 1)
        warning('Similarity values outside expected range [-1, 1]');
    end
end 