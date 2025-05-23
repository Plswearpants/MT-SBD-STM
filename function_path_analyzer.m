function function_path_analyzer()
    % Get the current script's path
    [currentScriptPath, ~, ~] = fileparts(mfilename('fullpath'));
    
    % Get all .m files in the current directory and subdirectories
    mFiles = dir(fullfile(currentScriptPath, '**', '*.m'));
    
    % Process each .m file
    for i = 1:length(mFiles)
        filePath = fullfile(mFiles(i).folder, mFiles(i).name);
        [~, fileName, ~] = fileparts(mFiles(i).name);
        
        % Read the file content with error handling
        try
            fid = fopen(filePath, 'r');
            if fid == -1
                fprintf('Warning: Could not open file %s\n', filePath);
                continue;
            end
            fileContent = fread(fid, '*char')';
            fclose(fid);
            
            % Find all function declarations
            functionMatches = regexp(fileContent, 'function\s+([a-zA-Z0-9_]+)', 'tokens');
            
            % Print the main function path
            fprintf('%s: %s\n', fileName, filePath);
            
            % Print subfunctions if any
            if ~isempty(functionMatches)
                for j = 2:length(functionMatches)  % Start from 2 to skip main function
                    subFuncName = functionMatches{j}{1};
                    fprintf('  └─ %s: %s\n', subFuncName, filePath);
                end
            end
        catch ME
            fprintf('Error processing file %s: %s\n', filePath, ME.message);
            if fid ~= -1
                fclose(fid);
            end
        end
    end
end 