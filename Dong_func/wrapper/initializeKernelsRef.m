function [data, params] = initializeKernelsRef(data, params, varargin)
%INITIALIZEKERNELSREF Wrapper for initializing kernels from reference slice
%
%   [data, params] = initializeKernelsRef(data, params, ...)
%
%   INPUTS:
%       data                - Data struct from previous block (e.g., from generateSyntheticData)
%       params              - Parameter struct from previous block
%
%       OPTIONAL (Name-Value pairs):
%       'kernel_selection_type' - 'selected' (interactive) or 'random' (default: 'selected')
%       'window_type'       - Window function to apply to kernels (default: {'gaussian', 2.5})
%                             Options: 'hann', 'hamming', 'blackman'
%                                      {'gaussian', alpha}
%                                      {'kaiser', beta}
%                                      {} or '' for no window
%
%   OUTPUTS:
%       data                - Updated data struct with new field:
%                             data.A1 - {1 × K} cell array of initialized kernels
%       params              - Updated parameter struct with initialization settings
%
%   DESCRIPTION:
%       This wrapper encapsulates the kernel initialization process for the
%       reference slice. It handles:
%       - Extracting kernel sizes for the reference slice
%       - Displaying ground truth kernels (only for 'selected' mode)
%       - Calling initialize_kernels() for interactive/random selection
%       - Displaying initialized kernels (always)
%       - Returning kernels with proper windowing applied
%
%   EXAMPLE:
%       % Interactive selection with Gaussian window
%       [data, params] = initializeKernelsRef(data, params, ...
%           'kernel_selection_type', 'selected', ...
%           'window_type', {'gaussian', 2.5});
%
%       % Random initialization without window
%       [data, params] = initializeKernelsRef(data, params, ...
%           'kernel_selection_type', 'random', ...
%           'window_type', {});
%
%   See also: initialize_kernels, apply_window, proj2oblique

    % Parse input arguments
    p = inputParser;
    addRequired(p, 'data', @isstruct);
    addRequired(p, 'params', @isstruct);
    addParameter(p, 'kernel_selection_type', 'selected', @(x) ismember(x, {'selected', 'random'}));
    addParameter(p, 'window_type', {'gaussian', 2.5});
    parse(p, data, params, varargin{:});
    
    % Extract parameters
    kernel_selection_type = p.Results.kernel_selection_type;
    window_type = p.Results.window_type;
    
    % Automatic display logic:
    % - Show ground truth ONLY for 'selected' mode (to guide user selection)
    % - Always show initialized kernels (to verify initialization)
    show_ground_truth = strcmp(kernel_selection_type, 'selected');
    show_initialized = true;
    
    % Validate required fields in data struct
    required_fields = {'Y_ref', 'A0_ref'};
    for i = 1:length(required_fields)
        if ~isfield(data, required_fields{i})
            error('Data struct must contain field: %s', required_fields{i});
        end
    end
    
    % Validate required fields in params struct
    required_params = {'num_kernels', 'ref_slice', 'kernel_sizes'};
    for i = 1:length(required_params)
        if ~isfield(params, required_params{i})
            error('Params struct must contain field: %s', required_params{i});
        end
    end
    
    % Extract kernel sizes for reference slice
    % params.kernel_sizes is [num_slices × num_kernels × 2]
    kernel_sizes_ref = reshape(params.kernel_sizes(params.ref_slice,:,:), [params.num_kernels, 2]);
    
    % Display ground truth kernels (if available and requested)
    if show_ground_truth && isfield(data, 'A0')
        fprintf('\n');
        fprintf('========================================================================\n');
        fprintf('  GROUND TRUTH(GT) KERNELS - Reference for Selection Order\n');
        fprintf('========================================================================\n');
        figure('Name', 'IN01A: Ground Truth Kernels (Reference Slice)');
        for k = 1:params.num_kernels
            subplot(1, params.num_kernels, k);
            imagesc(data.A0{params.ref_slice, k});
            axis square;
            title(sprintf('GT Kernel %d', k), 'FontWeight', 'bold', 'FontSize', 12);
            colorbar;
            colormap(gray);
        end
        sgtitle(sprintf('GT Kernels - Slice %d (SELECT IN THIS ORDER)', params.ref_slice), ...
                'FontWeight', 'bold', 'FontSize', 14);
        fprintf('  GT kernels are displayed above.\n');
        fprintf('========================================================================\n\n');
    end
    
    % Initialize kernels
    fprintf('  Initializing kernels for reference slice...\n');
    if strcmp(kernel_selection_type, 'selected')
        fprintf('\n');
        fprintf('  ╔═══════════════════════════════════════════════════════════════════╗\n');
        fprintf('  ║IMPORTANT: Select initialized kernels in the SAME ORDER as ground truth! ║\n');
        fprintf('  ║                                                                   ║\n');
        fprintf('  ║  → Look at the "Ground Truth Kernels" figure window above         ║\n');
        fprintf('  ║  → Select initialized kernel 1 first, then kernel 2, etc.         ║\n');
        fprintf('  ║  → This ensures correct metric calculation later                  ║\n');
        fprintf('  ╚═══════════════════════════════════════════════════════════════════╝\n\n');
    else
        fprintf('  Generating %d random initialized kernels.\n', params.num_kernels);
    end
    
    A1 = initialize_kernels(data.Y_ref, params.num_kernels, kernel_sizes_ref, ...
                            kernel_selection_type, window_type);
    
    % Display initialized kernels
    if show_initialized
        fprintf('  Displaying initialized kernels...\n');
        figure('Name', 'IN01A: Initialized Kernels (Reference Slice)');
        for n = 1:params.num_kernels
            subplot(1, params.num_kernels, n);
            imagesc(A1{n});
            title(sprintf('Initial Kernel %d', n));
            colorbar;
            axis square;
            colormap(gray);
        end
        sgtitle(sprintf('Initialized Kernels (Reference Slice %d)', params.ref_slice));
    end
    
    fprintf('  Kernel initialization complete.\n');
    
    % Add initialized kernels to data struct
    data.A1 = A1;
    
    % Update params struct with initialization settings
    params.kernel_selection_type = kernel_selection_type;
    params.window_type = window_type;
    params.kernel_sizes_ref = kernel_sizes_ref;
end

