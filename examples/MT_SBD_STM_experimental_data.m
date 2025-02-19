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

%% I.1Crop data
mask= gridMaskSquare(target_data);
target_data= gridCropMask(target_data, mask);
imagesc(target_data);
colormap("gray")
axis square

%% I.2 noise leveling 
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
A0 = randn([square_size n]);
% Need to put each slice back onto the sphere
A0 = proj2oblique(A0);

%% Initialize as specific crop of the target slice
[square_size,position, mask] = squareDrawSize(target_data(:,:,selected_slice));           	% determine kernel size
[A0, ~] = gridCropMask(target_data(:,:,selected_slice), mask);   % the cropped real data as kernel(wishlist, could act as initial kernel in the iteration process)
% Need to put each slice back onto the sphere
A0 = proj2oblique(A0);

%% Initialize as simulated kernel from TB model 
% Randomly choose n kernel slices from simulated LDoS data
n = 1;
load('example_data/LDoS_sim.mat');
sliceidx = randperm(size(LDoS_sim,3), n);
A0 = NaN([square_size n]);
    for i = 1:n
        A0 = imresize(LDoS_sim(:,:,sliceidx), square_size);
    end

% Need to put each slice back onto the sphere
A0 = proj2oblique(A0);
        
    
%% (Opt) Activation map generation:
% Generate activation map based on the sliced data
X0=activationCreateClick(target_data(:,:,selected_slice));

m = size(X0);          % image size for each slice / observation grid

%% (Opt)noise level determination 
eta_data = estimate_noise(target_data, 'std');  
SNR_data= var(A0(:))/eta_data;
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

% Create a VideoWriter object
%video_filename = 'Ag.avi';  % Specify the output file name
%v = VideoWriter(video_filename);  % Create the VideoWriter object
%v.FrameRate = 10;  % Set the frame rate (adjust as needed)
%open(v);  % Open the file for writing

% A function for showing updates as RTRM runs
figure;
dispfun = @( Y, A, X, square_size, kplus, idx ) showims_multikernel(Y,A0,X0,A,X,square_size,kplus,idx);

% Capture the current frame
frame = getframe(gcf);  % gcf gets the current figure
    
% Write the frame to the video
%writeVideo(v, frame);

% 2. The fun part
[Aout, Xout, extras] = SBD_test( Y, square_size, params, dispfun, A0);
%[Aout, Xout, extras] = SBD( Y, square_size, params, dispfun );

% close video
%close(v);

% Save the result
save('SBD-STM_datarun_Stephanie_AgGrid006_2.mat', 'Y', 'Xout', 'Aout','extras','square_size', "params");
%disp(['Video saved as ', video_filename]);
%% Visualization 
showims(Y,A0,Xout,Aout,Xout,square_size,[],1)

%% Visualization II
showims_fft(Y,A0,Xout,Aout,Xout,square_size,[],1)
