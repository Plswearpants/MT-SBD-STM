function [Aout, Xout, bout, extras, dataset_A0, dataset_X0] = adaptParallelResults(loaded_data)
    % ADAPTPARALLELRESULTS Converts parallel results format to visualization format
    %
    % Input:
    %   loaded_data: Structure containing the saved parallel results
    %
    % Outputs:
    %   Aout, Xout, bout, extras: Unwrapped data in 1xnum_kernel format
    %   dataset_A0, dataset_X0: Original dataset references
    
    % Unwrap the cell arrays
    Aout = loaded_data.Aout{1};
    Xout = loaded_data.Xout{1};
    bout = loaded_data.bout{1};
    extras = loaded_data.extras{1};
    dataset_A0 = loaded_data.dataset_A0{1};
    dataset_X0 = loaded_data.dataset_X0{1};
end
