function combined_score = computeCombined_activationScore(demixing_score, activation_recovery)
    % Compute combined score using geometric mean
    % Input:
    %   demixing_score: score from demixing metric [0,1]
    %   activation_recovery: score from activation recovery [0,1]
    % Output:
    %   combined_score: geometric mean of inputs [0,1]
    
    % Input validation
    if any(demixing_score(:) < 0 | demixing_score(:) > 1 | ...
           activation_recovery(:) < 0 | activation_recovery(:) > 1)
        error('Input metrics must be in range [0,1]');
    end
    
    % Compute geometric mean
    combined_score = sqrt(demixing_score .* activation_recovery);
end 