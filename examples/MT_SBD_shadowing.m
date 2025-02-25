% Description: This file is used to create a shadowing synthetic dataset
% according to an experimental dataset

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
[header, par, I, dIdV, LockindIdV, bias, midV, QPI, LockinQPI] = load3dsall('Grid Spectroscopy012.3ds', 5);
xsize = header.grid_dim(1);
ysize = header.grid_dim(2);
elayer = header.points;
estart = par(1);
eend = par(2);

%% ~~~~~~~~~~~~~~~~~~~~~~~~~Data preprocess~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
%% I.0 data selection 
target_data = LockindIdV;
rangetype='dynamic';
figure;
d3gridDisplay(target_data,rangetype);
selected_slice = input('Enter the slice number you want to analyze: ');
% change target data to single slice 
target_data = target_data(:,:,selected_slice);

%% I.1 (Opt) Crop data
mask= gridMaskSquare(target_data);
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
[square_size] = squareDrawSize(target_data);
%% Initialize as random kernel 
kerneltype = 'random';   
n = 1;               	% number of kernel slices
% Randomly generate n kernel slices
A1 = randn([square_size n]);
% Need to put each slice back onto the sphere
A1 = proj2oblique(A1);

%% Initialize as specific crop of the target slice
[square_size,position, mask] = squareDrawSize(target_data);           	% determine kernel size
[A1, ~] = gridCropMask(target_data, mask);   % the cropped real data as kernel(wishlist, could act as initial kernel in the iteration process)
% Need to put each slice back onto the sphere
A1 = proj2oblique(A1);

%% Initialize as simulated kernel from TB model 
% Randomly choose n kernel slices from simulated LDoS data
n = 1;
load('example_data/LDoS_sim.mat');
sliceidx = randperm(size(LDoS_sim,3), n);
A1 = NaN([square_size n]);
    for i = 1:n
        A1 = imresize(LDoS_sim(:,:,sliceidx), square_size);
    end

% Need to put each slice back onto the sphere
A1 = proj2oblique(A1);
        
    
%% (Opt) Activation map generation:
% Generate activation map based on the sliced data
X0=activationCreateClick(target_data(:,:,selected_slice));

m = size(X0);          % image size for each slice / observation grid

%% (Opt)noise level determination 
eta_data = estimate_noise(target_data, 'std');  
SNR_data= var(A1(:))/eta_data;
fprintf('SNR_data = %d', SNR_data);

%% (ESS) Define the observation as target_data
Y= target_data;

%% II. Sparse Blind Deconvolution:
%  ===============================
%% III parameter setting and SBD run and record the whole update into a video

% SBD settings
params.lambda1 = 0.1;              % regularization parameter for Phase I

params.phase2 = true;               % whether to do Phase II (refinement)
params.kplus = ceil(0.5 * square_size);       % padding for sphere lifting
params.lambda2 = 0.05;              % FINAL reg. param. value forclose  Phase II
params.nrefine = 3;                 % number of refinements

% Want entries of X to be nonnegative: see SBD_main.m
params.signflip = 0.2;
params.xpos     = true;
params.getbias  = true;
params.Xsolve = 'FISTA';

% A function for showing updates as RTRM runs
figure;
dispfun = @( Y, A, X, square_size, kplus, idx ) showims_multikernel(Y,A1,X,A,X,square_size,kplus,idx);


% 2. The fun part
[Aout, Xout, extras] = SBD_test( Y, square_size, params, dispfun, A1);

% Save the result
save('SBD-STM_datarun_Stephanie_AgGrid006_2.mat', 'Y', 'Xout', 'Aout','extras','square_size', "params");
%% Visualization 
showims(Y,A1,Xout,Aout,Xout,square_size,[],1)

%% Visualization II
showims_fft(Y,A1,Xout,Aout,Xout,square_size,[],1)
