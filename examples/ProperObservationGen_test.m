% Define required inputs 

SNR=2;
N_obs= 100;
observation_resolution=3;
defect_density=0.001;
LDoS_path = 'example_data/LDoS_sim.mat';

[Y, A0, X0, params] = properGenObservation(...
    SNR, ...
    N_obs, ...
    observation_resolution, ...
    defect_density);
    % Required Inputs:
    %   SNR: Signal-to-noise ratio
    %   N_obs: Size of observation lattice (N_obs x N_obs)
    %   observation_resolution: Resolution factor (pixels per lattice site)
    %   defect_density: Surface defect density (between 0 and 1)
    %
    % Optional Inputs:
    %   rho_single: Single defect QPI simulation at energy omega (3D array for multiple energies)
    %   LDoS_path: Path to LDoS simulation data

kernel_sizes = params.kernel_sizes;
A0_noiseless = params.A0_noiseless;
%% normalize and visualize
rangetype = 'dynamic';  
Y = normalizeBackgroundToZeroMean3D(Y,rangetype); 
Y = proj2oblique(Y);

figure; 
imagesc(Y);
colorbar;
axis square;

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
num_kernels = params.num_kernels;
kernel_sizes = params.kernel_sizes;
% Initialize kernels with the selected window type
A1 = initialize_kernels(Y, params.num_kernels, params.kernel_sizes, kerneltype, window_type);

% Display initialized kernels
figure;
for n = 1:params.num_kernels
    subplot(1, params.num_kernels, n);
    imagesc(A1{n});
    title(sprintf('Initial Kernel %d', n));
    colorbar;
    axis square;
end
%% 3. Settings

% A function for showing updates as RTRM runs
figure;

dispfun = cell(1,num_kernels);
dispfun{1} = @(Y, A, X, kernel_sizes, kplus) showims(Y,A0{1},X0(:,:,1),A,X,kernel_sizes,kplus,1); % here the last entry in the showims function is the energy layer index n. 
dispfun{2} = @(Y, A, X, kernel_sizes, kplus) showims(Y,A0{2},X0(:,:,2),A,X,kernel_sizes,kplus,1);
dispfun{3} = @(Y, A, X, kernel_sizes, kplus) showims(Y,A0{3},X0(:,:,3),A,X,kernel_sizes,kplus,1);

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

