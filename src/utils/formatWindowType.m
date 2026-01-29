function str = formatWindowType(window_type)
%FORMATWINDOWTYPE Convert window_type to string for logging
%   str = formatWindowType(window_type)
%
%   Handles various window_type formats:
%   - Cell array: {'gaussian', 2.5} -> '{gaussian, 2.50}'
%   - String/char: 'hann' -> 'hann'
%   - Empty: [] or {} -> 'none'
%
%   INPUTS:
%       window_type - Window type specification
%
%   OUTPUTS:
%       str - Formatted string representation

    if iscell(window_type) && ~isempty(window_type)
        if length(window_type) == 2
            str = sprintf('{%s, %.2f}', window_type{1}, window_type{2});
        else
            str = sprintf('{%s}', window_type{1});
        end
    elseif ischar(window_type) || isstring(window_type)
        str = char(window_type);
    elseif isempty(window_type)
        str = 'none';
    else
        str = 'custom';
    end
end

