function [Y, rho_single_resized_noisy, X_upsampled, params] = properGenObservation(varargin)
    % This function generates synthetic STM observations using a more physical approach
    % 
    % Required Inputs:
    %   SNR: Signal-to-noise ratio
    %   N_obs: Size of observation lattice (N_obs x N_obs)
    %   observation_resolution: Resolution factor (pixels per lattice site)
    %   defect_density: Surface defect density (between 0 and 1)
    %
    % Optional Inputs:
    %   rho_single: Single defect QPI simulation at energy omega (3D array for multiple energies)
    %   LDoS_path: Path to LDoS simulation data
    %
    % Outputs:
    %   Y: Final observation (pN_obs x pN_obs pixels)
    %   rho_single_resized: Cell array of resized single defect patterns
    %   X_upsampled: Upsampled activation maps (pN_obs x pN_obs x num_kernels)
    %   SNR: Signal-to-noise ratio used
    %   params: Structure containing all parameters used
    %
    % Example:
    %   % Basic usage - will prompt for LDoS file selection
    %   SNR = 5;              % Signal-to-noise ratio
    %   N_obs = 50;           % Observation lattice size (50x50)
    %   observation_resolution = 3;                % Resolution: 3 pixels per lattice site
    %   defect_density = 0.05;         % 5% defect density
    %   
    %   [Y, rho, X, SNR, params] = properGenObservation(SNR, N_obs, observation_resolution, defect_density);
    %   
    %   % The function will:
    %   % 1. Prompt for LDoS data file selection
    %   % 2. Display data for slice selection
    %   % 3. Allow user to select energy slices
    %   % 4. Process and display results
    
    % Parse inputs
    p = inputParser;
    p.addRequired('SNR');
    p.addRequired('N_obs');
    p.addRequired('observation_resolution');
    p.addRequired('defect_density');
    p.addParameter('rho_single', []);  % Optional: direct input of rho_single
    p.addParameter('LDoS_path', '');   % Optional: path to load data
    p.parse(varargin{:});
    
    % Extract parameters
    rho_single = p.Results.rho_single;
    SNR = p.Results.SNR;
    N_obs = p.Results.N_obs;
    p_scale = p.Results.observation_resolution;
    rho_d = p.Results.defect_density;
    LDoS_path = p.Results.LDoS_path;
    
    % Handle data loading
    if isempty(rho_single)
        if isempty(LDoS_path)
            [filename, pathname] = uigetfile({'*.mat'}, 'Select LDoS Result File');
            if isequal(filename, 0)
                error('No file selected');
            end
            LDoS_path = fullfile(pathname, filename);
        end
        
        try
            loaded_data = load(LDoS_path);
            % Load N_single from file
            if isfield(loaded_data, 'N')
                N_single = loaded_data.N;
            else
                error('N not found in the loaded file. This parameter is required.');
            end
            
            % Load rho_single
            if isfield(loaded_data, 'LDoS_sim')
                rho_single = loaded_data.LDoS_sim;
            else
                potential_fields = find(structfun(@(x) ndims(x) >= 2, loaded_data));
                field_names = fieldnames(loaded_data);
                
                if length(potential_fields) == 1
                    rho_single = loaded_data.(field_names{potential_fields});
                elseif length(potential_fields) > 1
                    [idx, ok] = listdlg('ListString', field_names(potential_fields), ...
                                      'SelectionMode', 'single', ...
                                      'PromptString', 'Select the QPI data field:');
                    if ~ok
                        error('No field selected');
                    end
                    rho_single = loaded_data.(field_names{potential_fields(idx)});
                else
                    error('No suitable data field found in the file');
                end
            end
        catch ME
            error('Error loading file: %s', ME.message);
        end
    else
        error('When providing rho_single directly, N_single must also be provided');
    end
    
    % Validate and prepare data
    if isempty(rho_single) || ~isnumeric(rho_single) || ndims(rho_single) < 2
        error('Invalid data format');
    end
    if ndims(rho_single) == 2
        rho_single = rho_single(:,:,1);
    end
    
    % Display full dataset for user review
    disp('Displaying full dataset. Use slider to review all slices.');
    disp('Close the display window when ready to proceed with slice selection.');
    figure;
    d3gridDisplay(rho_single, 'dynamic');
    
    % Handle slice selection
    total_slices = size(rho_single, 3);
    fprintf('\nTotal available slices: %d\n', total_slices);
    
    % Get slice selection from command line
    while true
        input_str = input('Enter slice numbers separated by spaces (e.g., "1 5 10"): ', 's');
        selected_indices = str2num(input_str); %#ok<ST2NM>
        
        % Validate input
        if isempty(selected_indices)
            fprintf('Invalid input. Please enter numbers separated by spaces.\n');
            continue;
        end
        
        if any(selected_indices < 1) || any(selected_indices > total_slices)
            fprintf('Invalid slice numbers. Must be between 1 and %d.\n', total_slices);
            continue;
        end
        
        if length(selected_indices) < 1
            fprintf('Please select at least one slice.\n');
            continue;
        end
        
        break;
    end
    
    % Sort indices and remove duplicates
    selected_indices = unique(sort(selected_indices));
    fprintf('Selected slices: %s\n', num2str(selected_indices));
    
    % Keep only selected slices
    rho_single = rho_single(:,:,selected_indices);
    selected_slices = selected_indices;
    num_kernels = length(selected_indices);
    
    % Initialize arrays for processed data
    rho_single_resized = cell(1, num_kernels);
    rho_single_resized_noisy = cell(1, num_kernels);
    X_upsampled = zeros(N_obs*p_scale, N_obs*p_scale, num_kernels);
    cutoff_M = zeros(1, num_kernels);
    target_size = zeros(num_kernels,2);
    % Process each selected slice
    for k = 1:num_kernels
        % Find cutoff M for this slice
        [M, M_pixels,~] = find_cutoff_noise_intersection(rho_single(:,:,k), SNR, N_single, false);
        cutoff_M(k) = M;
        
        % Get center and crop rho_single
        [ny, nx] = size(rho_single(:,:,k));
        center = ceil([ny, nx]/2);
        pixels_per_lattice = nx/N_single;
        
        % Crop the pattern
        range_y = max(1, center(1)-M_pixels):min(ny, center(1)+M_pixels);
        range_x = max(1, center(2)-M_pixels):min(nx, center(2)+M_pixels);
        rho_single_cutoff = rho_single(range_y, range_x, k);
        
        % Resize cutoff pattern to match observation scale
        target_size(k,:) = round([2*M*p_scale 2*M*p_scale]);
        rho_single_resized{k} = imresize(rho_single_cutoff, target_size(k,:));
        
        % Generate activation map with strict density control
        target_defects = round(N_obs * N_obs * rho_d);  % Expected number of defects
        tolerance = 0.1;  % 10% tolerance
        min_defects = max(round(target_defects * (1-tolerance)), 1);
        max_defects = round(target_defects * (1+tolerance));
        
        X_good = false;
        while ~X_good
            X = double(rand(N_obs, N_obs) <= rho_d);
            num_defects = sum(X(:));
            X_good = (num_defects >= min_defects) && (num_defects <= max_defects);
        end
        
        % Upsample X
        X_upsampled(:,:,k) = upsample_with_zero_blocks(X, p_scale);
    end
    
    % Display summary of processing
    figure('Name', 'Processing Summary');
    for k = 1:num_kernels
        subplot(3, num_kernels, k)
        imagesc(rho_single(:,:,k));
        title(sprintf('Original Slice %d', selected_slices(k)));
        axis square
        
        subplot(3, num_kernels, k + num_kernels)
        imagesc(rho_single_resized{k});
        title(sprintf('Processed (M=%.1f)', cutoff_M(k)));
        axis square
        
        subplot(3, num_kernels, k + 2*num_kernels)
        imagesc(X_upsampled(:,:,k));
        title('Activation Map');
        axis square
    end
    sgtitle('Processing Summary for Selected Slices');
    
    % Generate clean observation through convolution
    Y_clean = zeros(N_obs*p_scale, N_obs*p_scale);
    for k = 1:num_kernels
        Y_clean = Y_clean + convfft2(rho_single_resized{k}, X_upsampled(:,:,k));
    end
    
    % Add noise based on SNR
    eta = var(Y_clean(:)) / SNR;
    Y = Y_clean + sqrt(eta) * randn(size(Y_clean));

    for k =1:num_kernels
        rho_single_resized_noisy{k} = rho_single_resized{k} + sqrt(eta) * randn(size(rho_single_resized{k}));
    end

    % Package parameters
    params = struct();
    params.N_single = N_single;
    params.N_obs = N_obs;
    params.observation_resolution = p_scale;
    params.defect_density = rho_d;
    params.SNR = SNR;
    params.num_kernels = num_kernels;
    params.selected_slices = selected_slices;
    params.cutoff_M = cutoff_M;
    params.kernel_sizes = target_size;
    params.A0_noiseless = rho_single_resized;
    % Final visualization
    figure('Name', 'Final Observation');
    subplot(1,2,1)
    imagesc(Y_clean);
    colorbar;
    title('Clean Observation');
    axis square;
    
    subplot(1,2,2)
    imagesc(Y);
    colorbar;
    title(sprintf('Noisy Observation (SNR=%.1f)', SNR));
    axis square;
    
    sgtitle('Final Generated Observation');
end