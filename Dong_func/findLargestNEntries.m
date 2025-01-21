function [largest_values, row_indices, col_indices] = findLargestNEntries(A, n)
    % Flatten the matrix to a single vector
    A_flat = A(:);

    % Sort the vector in descending order and get the indices
    [sorted_A, sorted_indices] = sort(A_flat, 'descend');

    % Get the top n largest entries and their original indices
    largest_values = sorted_A(1:n);
    largest_indices = sorted_indices(1:n);

    % Convert linear indices to subscript indices
    [row_indices, col_indices] = ind2sub(size(A), largest_indices);
end
