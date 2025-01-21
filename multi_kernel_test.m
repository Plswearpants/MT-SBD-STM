clc; clear;

% Import Manopt and initialize the SBD package
run('../init_sbd');
fprintf('\n\n');

%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Create synthetic multi-kernel observation~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%% 0. Load the .3ds data
[header, par, I, dIdV, LockindIdV, bias, midV, QPI, LockinQPI] = load3dsall('Grid Spectroscopy006.3ds', 5);
xsize = header.grid_dim(1);
ysize = header.grid_dim(2);
elayer = header.points;
estart = par(1);
eend = par(2);

%% 0.1 Slice selection
d3gridDisplay(dIdV,'dynamic');
selected_slice = input('Enter the slice number you want to analyze: ');

%% I. SIMULATE DATA FOR SBD:
%  =========================

%% 1. Kernel settings
num_kernels = input('Enter the number of kernels to use: ');
n = 1;  % number of kernel slices per kernel

[square_size, position, mask] = squareDrawSize(dIdV(:,:,selected_slice));
[kernel_data, ~] = gridCropMask(dIdV(:,:,selected_slice), mask);

%% Initialize as simulated kernel from TB model 
% Randomly choose n kernel slices from simulated LDoS data
load('example_data/LDoS_sim.mat');
A0 = cell(1, num_kernels);
for k = 1:num_kernels
    sliceidx = randperm(size(LDoS_sim,3), n);
    A0{k} = NaN([square_size n]);
    for i = 1:n
        A0{k}(:,:,i) = imresize(LDoS_sim(:,:,sliceidx(i)), square_size);
    end
    % Need to put each slice back onto the sphere
    A0{k} = proj2oblique(A0{k});
end

%% 2. Activation map generation:
X0 = cell(1, num_kernels);
for k = 1:num_kernels
    X0{k} = activationCreateClick(dIdV(:,:,selected_slice));
end

m = size(X0{1});          % image size for each slice / observation grid

%% 3 normal noise level determination and generate eta_sim
eta_data = estimate_noise(dIdV(:,:,selected_slice),'std');  
SNR_data = var(kernel_data(:))/eta_data;
fprintf('SNR_data = %d\n', SNR_data);
SNR_sim = SNR_data;
eta_sim = var(A0{1}(:))/SNR_sim;

%% 4 observation generation:
Y = zeros([m n]);
for k = 1:num_kernels
    for i = 1:n
        Y(:,:,i) = Y(:,:,i) + convfft2(A0{k}(:,:,i), X0{k});     
    end
end
Y = Y + sqrt(eta_sim)*randn([m n]);
imagesc(Y);
colorbar;
axis square;

%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Start SBD-STM on the synthetic data~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%% Initialize as random kernel 
kerneltype = 'random';   

% Randomly generate n kernel slices
A1 = randn([square_size n]);
% Need to put each slice back onto the sphere
A1 = proj2oblique(A1);

%% 1. Settings

% A function for showing updates as RTRM runs
figure;
dispfun = @(Y, A, X, square_size, kplus, idx) showims(Y, A0, X0, A, X, square_size, kplus, idx);

% SBD settings
params.lambda1 = 1e-1;  % regularization parameter for Phase I

params.phase2 = true;
params.kplus = ceil(0.5 * square_size);
params.lambda2 = 5e-2;  % FINAL reg. param. value for Phase II
params.nrefine = 3;

params.signflip = 0.2;
params.xpos = true;
params.getbias = true;
params.Xsolve = 'FISTA';

% 2. The fun part
[Aout, Xout, extras] = SBD_test(Y, square_size, params, dispfun, A1);

% Save the result
save('SBD-STM.mat', 'Y', 'X0', 'A0', 'Xout', 'Aout', 'sliceidx', 'square_size');

%% Visualization I
figure();
showims_multi(Y, A0, X0, Aout, Xout, square_size, [], 5);

%% Visualization II
square_size = [size(Aout{1},1), size(Aout{1},1)];
showims_fft_multi(Y, A0, X0, Aout, Xout, square_size, [], 4, 1);
