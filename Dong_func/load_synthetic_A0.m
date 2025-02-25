function A0_all = load_synthetic_A0()
    % Initialize output cell array
    A0_all = cell(12, 1);  % 12 datasets
    
    % Define the path to the synthetic datasets
    folder_path = 'results/parallel results/full_synthetic_datasets_20241127_170302';
    
    % Load each dataset file
    for i = 1:12
        % Construct filename
        filename = sprintf('SBD_parallel_dataset%d_param_idx1.mat', i);
        filepath = fullfile(folder_path, filename);
        
        % Load data
        try
            loaded_data = load(filepath);
            % Extract A0 from the loaded data
            if isfield(loaded_data,'dataset_A0')
                A0_all{i} = loaded_data.dataset_A0{1,1};
            else
                warning('Could not find A0 in dataset %d', i);
                A0_all{i} = {};
            end
        catch ME
            warning('Error loading dataset %d: %s', i, ME.message);
            A0_all{i} = {};
        end
    end
    
    % Verify the structure
    fprintf('Loaded A0 from %d datasets:\n', length(A0_all));
    for i = 1:length(A0_all)
        if ~isempty(A0_all{i})
            fprintf('Dataset %d: A0 contains %d kernels\n', i, length(A0_all{i}));
        else
            fprintf('Dataset %d: No A0 found\n', i);
        end
    end
end 