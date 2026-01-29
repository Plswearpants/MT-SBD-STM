function scan_quality_factors(file_path)
    % Open the file
    fid = fopen(file_path, 'r');
    
    % Initialize variables
    quality_factors = [];
    iteration = 0;
    
    % Read the file line by line
    while ~feof(fid)
        line = fgetl(fid);
        
        % Use regular expression to find the quality factor
        [match, tokens] = regexp(line, 'Quality Factor = (\d+\.\d+e[-+]?\d+)', 'match', 'tokens');
        
        if ~isempty(match)
            iteration = iteration + 1;
            quality_factor = str2double(tokens{1}{1});
            quality_factors(iteration) = quality_factor;
            
            % Print each quality factor as it's found
            fprintf('Iteration %d: %e\n', iteration, quality_factor);
        end
    end
    
    % Close the file
    fclose(fid);
    
    % Print total number of quality factors found
    fprintf('\nTotal number of Quality Factors found: %d\n', length(quality_factors));
    
    % Plot the quality factors
    figure;
    semilogy(1:iteration, quality_factors, 'b-o');
    title('Quality Factor vs. Iteration');
    xlabel('Iteration');
    ylabel('Quality Factor');
    grid on;
end

% Usage example:
% scan_quality_factors('your_file.txt');