%% Block 1: Load the .3ds data

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
[header, par, I, dIdV, LockindIdV, bias, midV, QPI, LockinQPI] = load3dsall('Grid Spectroscopy002.3ds', 5);
xsize = header.grid_dim(1);
ysize = header.grid_dim(2);
elayer = header.points;
estart = par(1);
eend = par(2);
energy_range = linspace(estart, eend, elayer);
data_original = dIdV;
num_slices = size(data_original,3);
spatial = size(data_original,1);

%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Block 2: Data preprocessing
data_carried = data_original;
rangetype='dynamic';
figure;
d3gridDisplay(data_carried,rangetype);
slice_normalize = input('slice to normalize: ');
%% 2.2b defect masking
% Normalize background 
[data_carried] = normalizeBackgroundToZeroMean3D(data_carried, rangetype, slice_normalize);

f1=figure;
d3gridDisplay(data_carried,'dynamic');
index = input('Enter defect slice number: ');
num_defect_type = input('enter how many types of defects to mask: ');
close(f1);
% methods: 
% 1. Gaussian window "gw"
% 2. truncated gaussian gaussian smoothing "tg"
% 3. thresholding and remove defect features "threshold"
method = 'tg';

switch method
    case 'gw'
        % Apply Gaussian window masking
        [data_masked, ~] = defect_masking(data_carried, index);
    case 'tg'
        % Apply flat disk mask with Gaussian smoothing
        [data_masked, defect_loc] = gaussianMaskDefects(Y,index, num_defect_type);
    case 'threshold'
        % Apply threshold-based defect masking
        [data_masked, defect_mask] = thresholdDefects(data_carried, index);
    otherwise
        error('Unknown defect masking method. Choose "gw", "disk", or "threshold".');
end
data_carried = data_masked;

%% 2.3a: Correct streak
% Normalize background 
[data_carried] = normalizeBackgroundToZeroMean3D(data_carried, rangetype, slice_normalize);

[data_streakremoved, QPI_nostreaks] = RemoveStreaks(data_carried, 'Direction', 'vertical');
data_carried = data_streakremoved;

%% 2.4: crop dataset
mask= maskSquare(data_carried,0,slice_normalize);
data_cropped= gridCropMask(data_carried, mask);
data_carried = data_cropped;

%% 2.end: Normalize background 
[Y] = normalizeBackgroundToZeroMean3D(data_carried, rangetype, slice_normalize);

%% 3pre: Save the preprocessed data~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
save('Ag002_preprocessed_selected.mat', 'Y', 'data_masked', 'data_cropped', 'data_streakremoved', 'data_original', "energy_range",'energy_selected','defect_loc')


%% Block 3 data selection 
rangetype='dynamic';
f1=figure;
d3gridDisplay(Y,rangetype);
params.slices = input('input a list of slices: ');
Y=data_carried(:,:,params.slices);
energy_selected = energy_range(params.slices);
close(f1);
num_slices = size(Y,3);

%% 3.end: Normalize background 
[Y] = normalizeBackgroundToZeroMean3D(Y, rangetype, 1);

%% 4_pre Pick reference slice and Initialize reference kernels
f2=figure;
d3gridDisplay(Y,rangetype);
params.ref_slice = input('Enter reference slice number: ');
num_kernels = input('enter the number of kernels you wish to apply: ');
close(f2);

% Validate the input
if isempty(params.ref_slice) || ~isnumeric(params.ref_slice) || params.ref_slice < 1 || params.ref_slice > size(Y, 3)
    error('Invalid reference slice. Please enter a number between 1 and %d.', size(Y, 3));
end
fprintf('Using slice %d as reference\n', params.ref_slice);

Y = normalizeBackgroundToZeroMean3D(Y, 'dynamic', params.ref_slice);  % normalize Y
Y = proj2oblique(Y);

% Extract reference slice and initialize kernel
Y_ref = Y(:,:,params.ref_slice);
Y_ref = normalizeBackgroundToZeroMean3D(Y_ref, 'dynamic');  % normalize Y
Y_ref = proj2oblique(Y_ref);

% Display the selected reference slice
figure;
imagesc(Y_ref);
colorbar;
title(sprintf('Reference Slice %d', params.ref_slice));
axis square;

%% Initialize reference kernels
% draw square on the data to include as many visible ripples of the scattering as possible 
same_size = 1;
kerneltype = 'selected';
window_type = {'gaussian', 2.5};
%window_type = '';


if same_size
    [square_size] = squareDrawSize(Y_ref);
    %square_size=[80,80];
    kernel_sizes = repmat(square_size,[num_kernels,1]);
    A1_ref = initialize_kernels(Y_ref, num_kernels, kernel_sizes, kerneltype, window_type);
else
    A1_ref = cell(1, num_kernels);
    kernel_sizes = zeros(num_kernels, 2); % Store sizes of each kernel [height, width]
    for k = 1:num_kernels
        fprintf('Select region for kernel %d/%d\n', k, num_kernels);
        [square_size,position, mask] = squareDrawSize(Y_ref);           	% determine kernel size
        [A1_ref{k}, ~] = gridCropMask(Y_ref, mask);   % the cropped real data as kernel
        % Need to put each slice back onto the sphere
        A1_ref{k} = proj2oblique(A1_ref{k});
        % Store the kernel size
        kernel_sizes(k,:) = size(A1_ref{k});
    end
end

%% (ESS)noise level determination 
eta_data = estimate_noise(Y_ref, 'std');  

%% Block 4: Find Optimal Activation for Reference Slice
% Set up display functions
figure;
dispfun = cell(1,  num_kernels);
for n = 1:num_kernels
    dispfun{n} = @(Y, A, X, kernel_sizes, kplus) showims(Y_ref, A1_ref{n}, X, A, X, kernel_sizes, kplus, 1);
end

% SBD settings.
miniloop_iteration = 3;
outerloop_maxIT= 15;
%params_ref.energy = energy_selected(params.ref_slice);
params_ref.lambda1 = [0.05, 0.03,0.02, 0.05];  % regularization parameter for Phase I
%params_ref.lambda1 = [0.15, 0.15, 0.15, 0.15, 0.15];  % regularization parameter for Phase I
params_ref.phase2 = true;
params_ref.kplus = ceil(0.2 * kernel_sizes);
params_ref.lambda2 = [0.05, 0.05, 0.05, 0.15];  % FINAL reg. param. value for Phase II
params_ref.nrefine = 4;
params_ref.signflip = 0.2;
params_ref.xpos = true;
params_ref.getbias = true;
params_ref.Xsolve = 'FISTA';

% noise variance for computeResidualQuality.m
params_ref.noise_var = eta_data;

% Run and save 
% 2. The fun part
[A_ref, X_ref, b_ref, extras_ref] = MT_SBD(Y_ref, kernel_sizes, params_ref, dispfun, A1_ref, miniloop_iteration, outerloop_maxIT);

%% Visualize Reference result 
visualizeRealResult(Y_ref,A_ref, X_ref, b_ref, extras_ref);

%% Pad the A_ref to be size defined by user, normalize and use them as the A1
target_size = [100, 100];
kernel_sizes_pad = repmat(target_size,[num_kernels,1]);
%kernel_sizes_pad = [[120,120];[120,120];[65,65]];
A_pre_pad = A_ref;
A1_ref = cell(1, num_kernels);
for k = 1:num_kernels
    sz = size(A_pre_pad{k});
    pad_h = kernel_sizes_pad(k,1) - sz(1);
    pad_w = kernel_sizes_pad(k,2) - sz(2);

    % Calculate pre- and post-padding for centering
    pre_h = floor(pad_h / 2);
    post_h = ceil(pad_h / 2);
    pre_w = floor(pad_w / 2);
    post_w = ceil(pad_w / 2);

    % Pad so that the kernel is centered
    A1_ref{k} = padarray(A_pre_pad{k}, [pre_h, pre_w], 'pre');
    A1_ref{k} = padarray(A1_ref{k}, [post_h, post_w], 'post');
    
    A1_ref{k} = proj2oblique(A1_ref{k});
end

% visualize the padded kernels
figure;
for k = 1:num_kernels
    subplot(1,num_kernels,k);
    imagesc(A1_ref{k}); axis square;
    colorbar;
end

% define the intial activation map using X_ref
for k = 1:num_kernels
    params_ref.xinit{k}.X = X_ref(:,:,k);
    params_ref.xinit{k}.b = extras_ref.phase1.biter(k); 
end
%% Set ups before padded run 
% Set up display functions
figure;
dispfun = cell(1, num_kernels);
for n = 1:num_kernels
    dispfun{n} = @(Y, A, X, kernel_sizes, kplus) showims(Y_ref, A1_ref{n}, X, A, X, kernel_sizes, kplus, 1);
end

% SBD settings.
miniloop_iteration = 1;
outerloop_maxIT= 5;

params_ref.lambda1 = [0.03, 0.03, 0.03,0.03];  % regularization parameter for Phase I
%params_ref.lambda1 = [0.15, 0.15, 0.15, 0.15, 0.15];  % regularization parameter for Phase I
params_ref.phase2 = false;
params_ref.kplus = ceil(0.5 * kernel_sizes);
params_ref.lambda2 = [2e-2, 2e-2];  % FINAL reg. param. value for Phase II
params_ref.nrefine = 3;
params_ref.signflip = 0.2;
params_ref.xpos = true;
params_ref.getbias = true;
params_ref.Xsolve = 'FISTA';

% noise variance for computeResidualQuality.m
params_ref.noise_var = eta_data;
%% Run the padded initialization 
% 2. The fun part
[A_ref_pad, X_ref_pad, b_ref_pad, extras_ref_pad] = MT_SBD(Y_ref, kernel_sizes_pad, params_ref, dispfun, A1_ref, miniloop_iteration, outerloop_maxIT);

%% Visualize Padded result 
visualizeRealResult(Y_ref,A_ref_pad, X_ref_pad, b_ref_pad, extras_ref_pad);