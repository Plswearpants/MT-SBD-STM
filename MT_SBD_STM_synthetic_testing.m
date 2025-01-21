clc; clear;

% Import Manopt and initialize the SBD package
run('../init_sbd');
fprintf('\n\n');

%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Create synthetic multi-kernel observation~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%% 0. Define parameters
LDoS_path = 'example_data/LDoS_sim.mat';
interactive_selection = true;
num_kernels = 3;
image_size  = [300, 300];
kernel_sizes = zeros(num_kernels,2);
kernel_sizes(1,:) = [70, 70];
kernel_sizes(2,:) = [70, 70];
kernel_sizes(3,:) = [50, 50];
SNR = 5;
theta_cap = 2e-4;

%% 1. Generate synthetic observation and visualize
% generate 
[Y, A0_noiseless, A0, X0, SNR, params] = generateSyntheticSTMData(...
    num_kernels, ...
    image_size, ...
    kernel_sizes, ...
    SNR, ...
    theta_cap, ...
    interactive_selection, ...
    LDoS_path);

% reassign the parameters 
num_kernels = params.num_kernels;
kernel_sizes = params.kernel_sizes;
image_size = params.image_size;
SNR = params.SNR;

% normalize and visualize
rangetype = 'dynamic';  
Y = normalizeBackgroundToZeroMean3D(Y,rangetype); 
Y = proj2oblique(Y);

figure; 
imagesc(Y);
colorbar;
axis square;

%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Start SBD-STM on the synthetic data~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%% 2. Select starting kernel and display
kerneltype = 'selected';   % existing options: 'random' or 'selected'

% Initialize kernels with window option
%window_type = {};
 window_type = {'gaussian', 2};  % Example: gaussian window with alpha=2.5
% Other window options:
% window_type = 'hann';
% window_type = 'hamming';
% window_type = 'blackman';
% window_type = {'kaiser', 5};
% window_type = '';  % no window

% Initialize kernels with the selected window type
A1 = initialize_kernels(Y, params.num_kernels, params.kernel_sizes, kerneltype, window_type);

% Display initialized kernels
figure;
for n = 1:num_kernels
    subplot(1, num_kernels, n);
    imagesc(A1{n});
    title(sprintf('Initial Kernel %d', n));
    colorbar;
    axis square;
end
sgtitle('Initialized Kernels');
%% 3. Settings

% A function for showing updates as RTRM runs
figure;

dispfun = cell(1,num_kernels);
dispfun{1} = @(Y, A, X, kernel_size, kplus) showims(Y,A0{1},X0(:,:,1),A,X,kernel_size,kplus,1); % here the last entry in the showims function is the energy layer index n. 
dispfun{2} = @(Y, A, X, kernel_size, kplus) showims(Y,A0{2},X0(:,:,2),A,X,kernel_size,kplus,1);
dispfun{3} = @(Y, A, X, kernel_size, kplus) showims(Y,A0{3},X0(:,:,3),A,X,kernel_size,kplus,1);

% SBD settings.
initial_iteration = 3;
maxIT= 30;

params.lambda1 = [5e-2, 5e-2, 5e-2];  % regularization parameter for Phase I
params.phase2 = false;
params.kplus = ceil(0.5 * kernel_sizes);
params.lambda2 = [1e-2, 1e-2, 1e-2];  % FINAL reg. param. value for Phase II
params.nrefine = 3;
params.signflip = 0.2;
params.xpos = true;
params.getbias = true;
params.Xsolve = 'FISTA';

% this is for the test phase only:
params.X0 = X0;
params.A0 = A0;

%% Run and save 
% 2. The fun part
%[Aout, Xout, bout, extras] = SBD_test_multi(Y, kernel_size, params, dispfun, A1, initial_iteration, maxIT);
[Aout, Xout, bout, extras] = SBD_test_multi_demixing(Y, kernel_sizes, params, dispfun, A1, initial_iteration, maxIT);

% Save the result
% Generate a unique filename for the work space
timestamp = datestr(now, 'yyyymmdd_HHMMSS');
filename = sprintf('SBD_STM_results_%s.mat', timestamp);

% Ensure the filename is unique
counter = 1;
while exist(filename, 'file')
    counter = counter + 1;
    filename = sprintf('SBD_demixing_STM_results_3kernels%s_%d.mat', timestamp, counter);
end

% Save the specified variables to the workspace
save(filename, 'A0', 'X0', 'Aout', 'Xout', 'bout', 'extras', 'Y', 'SNR', 'A0_noiseless');

fprintf('Results saved to: %s\n', filename);

%% Visualization of Results
visualizeResults(Y, A0, Aout, X0, Xout, bout, extras);
%visualizeResults_old(Y, A0, Aout, X0, Xout, bout, extras);

