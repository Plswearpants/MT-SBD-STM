function update_config(config_file, param_path, new_value, new_file_name)
    % Check if config file exists
    if ~exist(config_file, 'file')
        error('Config file %s does not exist.', config_file);
    end
    
    % Load existing configuration
    config = load(config_file);
    
    % Update the specific parameter
    fields = strsplit(param_path, '.');
    current = config;
    for i = 1:length(fields)-1
        if ~isfield(current, fields{i})
            warning('Field %s not found. Creating it.', fields{i});
            current.(fields{i}) = struct();
        end
        current = current.(fields{i});
    end
    
    last_field = fields{end};
    if ~isfield(current, last_field)
        warning('Parameter %s not found in configuration. Adding it as a new parameter.', param_path);
    end
    current.(last_field) = new_value;
    
    % Update the main config structure
    if length(fields) > 1
        config = setfield(config, fields{:}, current.(last_field));
    else
        config.(fields{1}) = current.(last_field);
    end
    
    % Determine the file name to save to
    if nargin < 4 || isempty(new_file_name)
        save_file = config_file;
    else
        save_file = new_file_name;
    end
    
    % Save the updated configuration
    save(save_file, '-struct', 'config');
    
    if nargin >= 4 && ~isempty(new_file_name)
        fprintf('Updated configuration saved to %s\n', new_file_name);
    else
        fprintf('Existing configuration file updated: %s\n', config_file);
    end
end
