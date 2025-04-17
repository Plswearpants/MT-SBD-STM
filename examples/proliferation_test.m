%% Generate test set. 
SNR = 2;              % Signal-to-noise ratio
N_obs = 50;           % Observation lattice size (50x50)
observation_resolution = 4;  % Resolution: 3 pixels per lattice site
defect_density = 0.005;      % 0.1% defect density
%   
%   [Y, A0, X0, params] = properGen_full(SNR, N_obs, observation_resolution, defect_density);
%   
%   % Advanced usage with custom parameters
   [Y, A0, X0, params] = properGen_full(SNR, N_obs, observation_resolution, defect_density, ...
       'LDoS_path', 'example_data/LDoS_single_defect_save.mat', ...
       'num_slices', 15);

%% 