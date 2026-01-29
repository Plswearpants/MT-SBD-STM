function [best_slice_idx, similarity_scores] = findMostSimilarSlice(reference, volume, show_dist)
    % Find the most similar slice in a 3D volume to a 2D reference image
    % Automatically resizes volume slices to match reference dimensions
    %
    % Inputs:
    %   reference: 2D array (M×N)
    %   volume: 3D array (P×Q×K) where K is number of slices
    %   show_dist: boolean to show distribution plot (default: true)
    %
    % Outputs:
    %   best_slice_idx: Index of most similar slice
    %   similarity_scores: Array of similarity scores for all slices
    
    if nargin < 3
        show_dist = true;
    end
    
    % Input validation
    validateattributes(reference, {'numeric'}, {'2d'});
    validateattributes(volume, {'numeric'}, {'3d'});
    
    % Get dimensions
    ref_size = size(reference);
    vol_size = size(volume);
    num_slices = vol_size(3);
    
    % Resize if needed
    if ~isequal(ref_size, vol_size(1:2))
        volume_resized = zeros(ref_size(1), ref_size(2), num_slices);
        for i = 1:num_slices
            volume_resized(:,:,i) = imresize(volume(:,:,i), ref_size);
        end
    else
        volume_resized = volume;
    end
    
    % Normalize reference
    reference_norm = reference / norm(reference(:));
    
    % Compare each slice
    similarity_scores = zeros(num_slices, 1);
    for i = 1:num_slices
        slice = volume_resized(:,:,i);
        slice_norm = slice / norm(slice(:));
        similarity_scores(i) = abs(sum(reference_norm(:) .* slice_norm(:)));
    end
    
    % Find best match
    [max_score, best_slice_idx] = max(similarity_scores);
    
    % Plot distribution if requested
    if show_dist
        figure('Name', 'Similarity Score Distribution');
        
        % Histogram
        subplot(2,1,1);
        histogram(similarity_scores, min(20,num_slices), 'Normalization', 'probability');
        hold on;
        plot([max_score max_score], ylim, 'r--', 'LineWidth', 2);
        title('Distribution of Similarity Scores');
        xlabel('Similarity Score');
        ylabel('Probability');
        legend('Score Distribution', 'Best Match');
        
        % Bar plot
        subplot(2,1,2);
        bar(similarity_scores);
        hold on;
        plot([best_slice_idx best_slice_idx], ylim, 'r--', 'LineWidth', 2);
        title('Similarity Scores by Slice');
        xlabel('Slice Index');
        ylabel('Similarity Score');
        legend('Scores', 'Best Match');
        
        % Add text information
        text(0.02, 0.98, sprintf('Best Match: Slice %d\nScore: %.4f', ...
             best_slice_idx, max_score), ...
             'Units', 'normalized', ...
             'VerticalAlignment', 'top', ...
             'BackgroundColor', 'w');
    end
end 