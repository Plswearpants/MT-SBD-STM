function [param_sets, param_descriptions] = sample_parameter_space_general(param_ranges, varargin)
    % Sample parameter space for synthetic dataset testing
    %
    % Inputs:
    %   param_ranges: struct with fields containing parameter ranges
    %       .theta_cap: array of theta_cap values
    %       .kernel_size: array of kernel size values
    %       .SNR: array of SNR values
    %
    % Optional Input:
    %   sampling_style: string, default='full'
    %       'full': Uses all combinations of parameters
    %       'custom': For future custom sampling strategies
    %
    % Outputs:
    %   param_sets: [N x 3] array where each row is [theta_cap, kernel_size, SNR]
    %   param_descriptions: cell array of strings describing each parameter set
    
    % Input parsing
    p = inputParser;
    addRequired(p, 'param_ranges', @isstruct);
    addOptional(p, 'sampling_style', 'full', @ischar);
    parse(p, param_ranges, varargin{:});
    
    sampling_style = p.Results.sampling_style;
    
    % Validate param_ranges structure
    validateParamRanges(param_ranges);
    
    % Select sampling method
    switch sampling_style
        case 'full'
            [param_sets, param_descriptions] = full_parameter_sampling(param_ranges);
        case 'custom'
            % Add custom sampling methods here if needed
            error('Custom sampling not yet implemented');
        otherwise
            error('Unsupported sampling style: %s', sampling_style);
    end
end

function [param_sets, param_descriptions] = full_parameter_sampling(D)
    % Create all combinations of parameters using meshgrid
    [T, K, S] = meshgrid(D.theta_cap, D.kernel_size, D.SNR);
    
    % Reshape to get parameter sets
    param_sets = [T(:), K(:), S(:)];
    
    % Create descriptions
    num_sets = size(param_sets, 1);
    param_descriptions = cell(num_sets, 1);
    
    for i = 1:num_sets
        param_descriptions{i} = sprintf('θ=%.2e, A_ratio=%.3f, SNR=%.1f', ...
            param_sets(i,1), param_sets(i,2), param_sets(i,3));
    end
end

function validateParamRanges(param_ranges)
    required_fields = {'theta_cap', 'kernel_size', 'SNR'};
    
    for i = 1:length(required_fields)
        field = required_fields{i};
        assert(isfield(param_ranges, field), ...
            sprintf('param_ranges missing required field: %s', field));
        assert(~isempty(param_ranges.(field)), ...
            sprintf('%s must not be empty', field));
        assert(isnumeric(param_ranges.(field)), ...
            sprintf('%s must be numeric', field));
    end
end 