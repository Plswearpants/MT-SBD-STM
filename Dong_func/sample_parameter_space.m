function [param_sets, param_descriptions] = sample_parameter_space(param_ranges, varargin)
    % Sample parameter space for synthetic dataset testing using specified strategy
    %
    % Inputs:
    %   param_ranges: struct with fields containing parameter ranges [hard, medium, easy]
    %       .theta_cap
    %       .kernel_size
    %       .SNR
    %
    % Optional Input:
    %   sampling_style: string, default='targeted_sparse'
    %       'targeted_sparse': Focuses on edge cases and medium combinations
    %       (other styles could be added: 'full_factorial', 'latin_hypercube', etc.)
    %
    % Outputs:
    %   param_sets: [N x 3] array where each row is [theta_cap, kernel_size, SNR]
    %   param_descriptions: cell array of strings describing each parameter set
    
    % Input parsing
    p = inputParser;
    addRequired(p, 'param_ranges', @isstruct);
    addOptional(p, 'sampling_style', 'targeted_sparse', @ischar);
    parse(p, param_ranges, varargin{:});
    
    sampling_style = p.Results.sampling_style;
    
    % Validate param_ranges structure
    validateParamRanges(param_ranges);
    
    % Select sampling method
    switch sampling_style
        case 'targeted_sparse'
            [param_sets, param_descriptions] = targeted_sparse_sampling(param_ranges);
        % Add more sampling methods here as needed:
        % case 'full_factorial'
        %     [param_sets, param_descriptions] = full_factorial_sampling(param_ranges);
        % case 'latin_hypercube'
        %     [param_sets, param_descriptions] = latin_hypercube_sampling(param_ranges);
        otherwise
            error('Unsupported sampling style: %s', sampling_style);
    end
end

function [param_sets, param_descriptions] = targeted_sparse_sampling(D)
    % Initialize storage
    param_sets = zeros(12, 3);  % 12 sets × 3 parameters
    param_descriptions = cell(12, 1);
    set_idx = 1;
    
    % Group 1: Hardest case (all hard)
    param_sets(set_idx,:) = [D.theta_cap(1), D.kernel_size(1), D.SNR(1)];
    param_descriptions{set_idx} = 'Hardest case: dense activation, large kernel, low SNR';
    set_idx = set_idx + 1;
    
    % Group 2: Hard single parameter cases
    % Dense activation, small kernel, high SNR
    param_sets(set_idx,:) = [D.theta_cap(1), D.kernel_size(3), D.SNR(3)];
    param_descriptions{set_idx} = 'Hard theta: dense activation, small kernel, high SNR';
    set_idx = set_idx + 1;
    
    % Sparse activation, large kernel, high SNR
    param_sets(set_idx,:) = [D.theta_cap(3), D.kernel_size(1), D.SNR(3)];
    param_descriptions{set_idx} = 'Hard kernel: sparse activation, large kernel, high SNR';
    set_idx = set_idx + 1;
    
    % Sparse activation, small kernel, low SNR
    param_sets(set_idx,:) = [D.theta_cap(3), D.kernel_size(3), D.SNR(1)];
    param_descriptions{set_idx} = 'Hard SNR: sparse activation, small kernel, low SNR';
    set_idx = set_idx + 1;
    
    % Group 3: Easiest case (all easy)
    param_sets(set_idx,:) = [D.theta_cap(3), D.kernel_size(3), D.SNR(3)];
    param_descriptions{set_idx} = 'Easiest case: sparse activation, small kernel, high SNR';
    set_idx = set_idx + 1;
    
    % Group 4: Medium cases
    % All medium
    param_sets(set_idx,:) = [D.theta_cap(2), D.kernel_size(2), D.SNR(2)];
    param_descriptions{set_idx} = 'Medium all: medium activation, medium kernel, medium SNR';
    set_idx = set_idx + 1;
    
    % Medium-Medium-Easy combinations
    param_sets(set_idx,:) = [D.theta_cap(2), D.kernel_size(2), D.SNR(3)];
    param_descriptions{set_idx} = 'Medium-Medium-Easy: medium activation, medium kernel, high SNR';
    set_idx = set_idx + 1;
    
    % Medium-Easy-Medium combinations
    param_sets(set_idx,:) = [D.theta_cap(2), D.kernel_size(3), D.SNR(2)];
    param_descriptions{set_idx} = 'Medium-Easy-Medium: medium activation, small kernel, medium SNR';
    set_idx = set_idx + 1;
    
    % Medium-Easy-Easy combinations
    param_sets(set_idx,:) = [D.theta_cap(2), D.kernel_size(3), D.SNR(3)];
    param_descriptions{set_idx} = 'Medium-Easy-Easy: medium activation, small kernel, high SNR';
    set_idx = set_idx + 1;
    
    % Easy-Medium-Medium combinations
    param_sets(set_idx,:) = [D.theta_cap(3), D.kernel_size(2), D.SNR(2)];
    param_descriptions{set_idx} = 'Easy-Medium-Medium: sparse activation, medium kernel, medium SNR';
    set_idx = set_idx + 1;
    
    % Easy-Medium-Easy combinations
    param_sets(set_idx,:) = [D.theta_cap(3), D.kernel_size(2), D.SNR(3)];
    param_descriptions{set_idx} = 'Easy-Medium-Easy: sparse activation, medium kernel, high SNR';
    set_idx = set_idx + 1;
    
    % Easy-Easy-Medium combinations
    param_sets(set_idx,:) = [D.theta_cap(3), D.kernel_size(3), D.SNR(2)];
    param_descriptions{set_idx} = 'Easy-Easy-Medium: sparse activation, small kernel, medium SNR';
end

function validateParamRanges(param_ranges)
    required_fields = {'theta_cap', 'kernel_size', 'SNR'};
    
    for i = 1:length(required_fields)
        field = required_fields{i};
        assert(isfield(param_ranges, field), ...
            sprintf('param_ranges missing required field: %s', field));
        assert(length(param_ranges.(field)) == 3, ...
            sprintf('%s must have exactly 3 values [hard, medium, easy]', field));
    end
end 