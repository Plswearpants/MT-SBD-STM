function B = old2new_metric(A)
    % Convert old metric format (cell array of structs) to new format (n×2 array)
    %
    % Input:
    %   A: 1×n cell array, where each entry is a struct containing:
    %      - similarity: 1×2 array of similarity scores for each kernel
    %      (and potentially other fields that are not used)
    %
    % Output:
    %   B: n×2 array where:
    %      - Each row represents one iteration
    %      - Each column represents similarity scores for one kernel
    %
    % Example:
    %   A = {struct('similarity',[0.8 0.7]), struct('similarity',[0.85 0.75])}
    %   B = old2new_metric(A)  % Returns: [0.8 0.7; 0.85 0.75]
    
    n = length(A);      % number of iterations
    B = zeros(n, 2);    % preallocate output array for 2 kernels
    
    % Extract similarity scores from each iteration
    for i = 1:n
        B(i,:) = A{i}.similarity;  % get similarity field from struct
    end
end 