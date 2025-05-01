% This script runs MT_SBD on real dataset. 
clc; clear;

% Import Manopt and initialize the SBD package
run('../init_sbd');
fprintf('\n\n');

%% 0. Load the .3ds data

% INPUTS
% 1: Data file to load, including file type ('QPI.3ds' for example)
% 2: Smoothing sigma for current data

% OUTPUTS
% header: Variable containing all experimental parameters
% I: Current data, smoothed by sigma
% dIdV: Numerically differentiated current data
% voltage: Vector of voltages for current
% midV: Vector on voltages for dIdV/QPI (midpoint of voltage vector)
% QPI: Fourier transformed dIdV data

% Modified function load3dsall from supplied matlab code from Nanonis
[header, par, I, dIdV, LockindIdV, bias, midV, QPI, LockinQPI] = load3dsall('Grid Spectroscopy_QPI_24nm_3days002.3ds', 5);
xsize = header.grid_dim(1);
ysize = header.grid_dim(2);
elayer = header.points;
estart = par(1);
eend = par(2);

%% ~~~~~~~~~~~~~~~~~~~~~~~~~Data preprocess~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
%% I.0 data selection 
target_data = dIdV;
rangetype='dynamic';
figure;
d3gridDisplay(target_data,rangetype);
selected_slice = input('Enter the slice number you want to analyze: ');
num_kernels = input('enter the number of kernels you wish to apply: ');
% change target data to single slice 
target_data = target_data(:,:,selected_slice);

%% direct dataselect
target_data = data_filtered;
num_kernels=2;
rangetype='dynamic';
%% I.1 (Opt) Crop data
mask= maskRectangle(target_data);
target_data= gridCropMask(target_data, mask);
imagesc(target_data);
colormap("gray")
axis square

%% I.2 (Opt) noise leveling(if you have uneven structured noise)
% noise level determination 
eta_data = estimate_noise(target_data,'std');  

% level noise at a sepecific dimension to remove the streak noise 
target_data = levelNoiseInteractive(target_data,'x');

%% I.3 data normalization
target_data = normalizeBackgroundToZeroMean3D(target_data,rangetype); 

target_data = proj2oblique(target_data);

figure; 
imagesc(target_data);
colorbar;
axis square;
%% I. SIMULATE DATA FOR SBD:
%  =========================

%% ~~~~~~~~~~~Initialize Kernel guess - see kernel types below~~~~~~~~~~~~~

%% define square_size
% draw square on the data to include as many visible ripples of the scattering as possible 
same_size = 1;
kerneltype = 'random';
window_type = {'gaussian', 3};
%window_type = '';

if same_size
    [square_size] = squareDrawSize(target_data);
    kernel_sizes = repmat(square_size,[num_kernels,1]);
    A1 = initialize_kernels(target_data, num_kernels, kernel_sizes, kerneltype, window_type);
else
    A1 = cell(1, num_kernels);
    kernel_sizes = zeros(num_kernels, 2); % Store sizes of each kernel [height, width]
    for k = 1:num_kernels
        fprintf('Select region for kernel %d/%d\n', k, num_kernels);
        [square_size,position, mask] = squareDrawSize(target_data);           	% determine kernel size
        [A1{k}, ~] = gridCropMask(target_data, mask);   % the cropped real data as kernel
        % Need to put each slice back onto the sphere
        %A1{k} = proj2oblique(A1{k});
        % Store the kernel size
        kernel_sizes(k,:) = size(A1{k});
    end
end

%% Display initialized kernels
figure;
for n = 1:num_kernels
    subplot(1, num_kernels, n);
    imagesc(A1{n});
    title(sprintf('Initial Kernel %d', n));
    colorbar;
    axis square;
end
sgtitle('Initialized Kernels');
%% (ESS)noise level determination 
eta_data = estimate_noise(target_data, 'std');  

%% (Opt) determine SNR
SNR_data= var(A1{1,1})/eta_data;
fprintf('SNR_data = %d', SNR_data);

%% (ESS) Define the observation as target_data
Y= target_data;

%% II. Sparse Blind Deconvolution:
%  ===============================
%% Settings

% A function for showing updates as RTRM runs
figure;

dispfun = cell(1,num_kernels);
dispfun{1} = @(Y, A, X, kernel_sizes, kplus) showims(Y,A1{1},X,A,X,kernel_sizes,kplus,1); % here the last entry in the showims function is the energy layer index n. 
dispfun{2} = @(Y, A, X, kernel_sizes, kplus) showims(Y,A1{2},X,A,X,kernel_sizes,kplus,1);
%dispfun{3} = @(Y, A, X, kernel_sizes, kplus) showims(Y,A1{3},X,A,X,kernel_sizes,kplus,1);
%dispfun{4} = @(Y, A, X, kernel_sizes, kplus) showims(Y,A1{4},X,A,X,kernel_sizes,kplus,1);
%dispfun{5} = @(Y, A, X, kernel_sizes, kplus) showims(Y,A1{5},X,A,X,kernel_sizes,kplus,1);


% SBD settings.
miniloop_iteration = 2;
outerloop_maxIT= 30;

params.lambda1 = [0.1, 0.1];  % regularization parameter for Phase I
%params.lambda1 = [0.15, 0.15, 0.15, 0.15, 0.15];  % regularization parameter for Phase I
params.phase2 = false;
params.kplus = ceil(0.5 * kernel_sizes);
params.lambda2 = [2e-2, 2e-2];  % FINAL reg. param. value for Phase II
params.nrefine = 3;
params.signflip = 0.2;
params.xpos = true;
params.getbias = true;
params.Xsolve = 'FISTA';

% noise variance for computeResidualQuality.m
params.noise_var = eta_data;

%% Run and save 
% 2. The fun part
[Aout, Xout, bout, extras] = MT_SBD(Y, kernel_sizes, params, dispfun, A1, miniloop_iteration, outerloop_maxIT);

% Save the result
% Generate a unique filename for the work space
timestamp = datestr(now, 'yyyymmdd_HHMMSS');
filename = sprintf('grid012_SBD_STM_results_%s.mat', timestamp);

% Ensure the filename is unique
counter = 1;
while exist(filename, 'file')
    counter = counter + 1;
    filename = sprintf('SBD_demixing_STM_results_3kernels%s_%d.mat', timestamp, counter);
end

% Save the specified variables to the workspace
save(filename, 'Y', 'A1', 'Aout', 'Xout', 'bout', 'extras', 'Y');

fprintf('Results saved to: %s\n', filename);

%% Visualization of Results
visualizeResults(Y, A1, Aout, Xout, Xout, bout, extras);
%visualizeResults_old(Y, A0, Aout, X0, Xout, bout, extras);
