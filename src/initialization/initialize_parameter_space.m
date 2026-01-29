function params = initialize_parameter_space(dataset_ranges, system_ranges, fixed_params)
    % Initialize parameter space for synthetic dataset testing
    %
    % Inputs:
    %   dataset_ranges: struct with fields theta_cap, kernel_size, SNR
    %                  each field should contain a 3-element vector [hard, medium, easy]
    %   system_ranges: struct with fields mini_loop, lambda1
    %                 each field should contain a 3-element vector
    %   fixed_params: struct with fields for fixed parameters
    %
    % Example usage:
    %   dataset_ranges.theta_cap = [1e-3, 1e-4, 1e-5];
    %   dataset_ranges.kernel_size = [0.16, 0.04, 0.01];
    %   dataset_ranges.SNR = [1.0, 3.16, 10.0];
    %   system_ranges.mini_loop = [1, 3, 9];
    %   system_ranges.lambda1 = [1e-2, 3.16e-2, 1e-1];
    %   fixed_params.num_kernels = 2;
    %   fixed_params.n = 1;
    %   fixed_params.relative_theta_cap = 1;
    %   fixed_params.relative_kernel_size = 1;
    %   fixed_params.image_size = [300, 300];
    
    % Validate inputs
    validateInputs(dataset_ranges, system_ranges, fixed_params);
    
    % Store the parameter ranges
    params.dataset_params = dataset_ranges;
    params.system_params = system_ranges;
    params.fixed_params = fixed_params;
    
    % Generate test sets
    test_sets = generate_test_sets(dataset_ranges);
    params.test_sets = test_sets;
    
    % Calculate total number of test combinations
    params.num_test_sets = length(test_sets);
    params.num_system_combinations = length(system_ranges.mini_loop) * length(system_ranges.lambda1);
    params.total_tests = params.num_test_sets * params.num_system_combinations;
    
    % Add metadata
    params.metadata.creation_date = datetime('now');
    params.metadata.description = 'Parameter space for synthetic dataset testing';
end

function test_sets = generate_test_sets(D)
    % Initialize test sets structure
    test_sets = struct();
    set_idx = 1;
    
    % Group 1: Hardest case
    test_sets(set_idx) = create_test_set([D.theta_cap(1), D.kernel_size(1), D.SNR(1)], ...
        'hard-hard-hard', 'dense activation, large kernel, low SNR');
    set_idx = set_idx + 1;
    
    % Group 2: Hard single parameter cases
    hard_easy_combinations = [
        1, 3, 3;  % hard theta_cap, easy kernel_size, easy SNR
        3, 1, 3;  % easy theta_cap, hard kernel_size, easy SNR
        3, 3, 1   % easy theta_cap, easy kernel_size, hard SNR
    ];
    
    for i = 1:size(hard_easy_combinations, 1)
        params = [
            D.theta_cap(hard_easy_combinations(i,1)), ...
            D.kernel_size(hard_easy_combinations(i,2)), ...
            D.SNR(hard_easy_combinations(i,3))
        ];
        test_sets(set_idx) = create_test_set(params, ...
            sprintf('case-%d', set_idx), ...
            sprintf('Test set %d', set_idx));
        set_idx = set_idx + 1;
    end
    
    % Group 3: Easiest case
    test_sets(set_idx) = create_test_set([D.theta_cap(3), D.kernel_size(3), D.SNR(3)], ...
        'easy-easy-easy', 'sparse activation, small kernel, high SNR');
    set_idx = set_idx + 1;
    
    % Group 4: Medium cases
    % Generate all combinations of medium and easy (excluding all-easy)
    for i = 2:3  % theta_cap
        for j = 2:3  % kernel_size
            for k = 2:3  % SNR
                if ~(i == 3 && j == 3 && k == 3)  % Skip all-easy case
                    params = [D.theta_cap(i), D.kernel_size(j), D.SNR(k)];
                    test_sets(set_idx) = create_test_set(params, ...
                        sprintf('case-%d', set_idx), ...
                        sprintf('Test set %d', set_idx));
                    set_idx = set_idx + 1;
                end
            end
        end
    end
end

function test_set = create_test_set(params, difficulty, description)
    test_set = struct();
    test_set.params = params;
    test_set.difficulty = difficulty;
    test_set.description = description;
end

function validateInputs(dataset_ranges, system_ranges, fixed_params)
    % Validate that all required fields exist and have correct dimensions
    required_dataset_fields = {'theta_cap', 'kernel_size', 'SNR'};
    required_system_fields = {'mini_loop', 'lambda1'};
    required_fixed_fields = {'num_kernels', 'n', 'relative_theta_cap', ...
        'relative_kernel_size', 'image_size'};
    
    % Check dataset ranges
    for i = 1:length(required_dataset_fields)
        field = required_dataset_fields{i};
        assert(isfield(dataset_ranges, field), ...
            sprintf('dataset_ranges missing required field: %s', field));
        assert(length(dataset_ranges.(field)) == 3, ...
            sprintf('%s must have exactly 3 values', field));
    end
    
    % Check system ranges
    for i = 1:length(required_system_fields)
        field = required_system_fields{i};
        assert(isfield(system_ranges, field), ...
            sprintf('system_ranges missing required field: %s', field));
        assert(length(system_ranges.(field)) == 3, ...
            sprintf('%s must have exactly 3 values', field));
    end
    
    % Check fixed parameters
    for i = 1:length(required_fixed_fields)
        field = required_fixed_fields{i};
        assert(isfield(fixed_params, field), ...
            sprintf('fixed_params missing required field: %s', field));
    end
end

% Helper function to generate parameter combinations
function combinations = generate_parameter_combinations(test_set, system_params)
    [ML, L1] = ndgrid(system_params.mini_loop, system_params.lambda1);
    combinations = [ML(:), L1(:)];
end 