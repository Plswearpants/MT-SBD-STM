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
[header, par, I, dIdV, LockindIdV, bias, midV, QPI, LockinQPI] = load3dsall('QPImap012.3ds', 10);
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
% initialize the preprocessing parameters
preprocessing_params = struct();
data_carried = data_original;
rangetype='dynamic';
figure;
d3gridDisplay(data_carried,rangetype);
preprocessing_params.slice_normalize = input('slice to normalize: ');

%% 2.1: Remove bragg peaks
% Normalize background 
[data_carried] = normalizeBackgroundToZeroMean3D(data_carried, rangetype, preprocessing_params.slice_normalize);
% Bragg remove
[data_braggremoved]=removeBragg(data_carried);
data_carried = data_braggremoved;

%% 2.2: crop dataset
mask= maskSquare(data_carried,0,40,'square');
data_cropped= gridCropMask(data_carried, mask);
data_carried = data_cropped;

% Use full energy range; no slice selection needed for local streak removal
energy_selected = energy_range;

%% 2.3: Local streak removal using single reference slice workflow
interactive = true;
opts = struct();
[data_carried, preprocessing_params.streak_params] = streakRemovalWorkflow(data_carried, interactive, opts);

%% 2.4: Interpolation using single reference slice workflow
interactive_interp = true;
opts_interp = struct();
[data_carried, preprocessing_params.interp_params] = interpRemovalWorkflow(data_carried, interactive_interp, opts_interp);

%% 2.5 defect masking
% Normalize background 
[data_carried] = normalizeBackgroundToZeroMean3D(data_carried, rangetype, preprocessing_params.slice_normalize);

f1=figure;
d3gridDisplay(data_carried,'dynamic');
preprocessing_params.defect_slice = input('Enter defect slice number: ');
preprocessing_params.num_defect_type = input('enter how many types of defects to mask: ');
close(f1);
% methods: 
% 1. Gaussian window "gw"
% 2. truncated gaussian gaussian smoothing "tg"
% 3. thresholding and remove defect features "threshold"
preprocessing_params.defect_masking_method = 'tg';

switch preprocessing_params.defect_masking_method
    case 'gw'
        % Apply Gaussian window masking
        [data_masked, ~] = defect_masking(data_carried, preprocessing_params.defect_slice);
    case 'tg'
        % Apply flat disk mask with Gaussian smoothing
        % Interactive mask creation and application:
        %if isfield(preprocessing_params, 'defect_mask') && ~isempty(preprocessing_params.defect_mask)
            %[data_masked, ~] = gaussianMaskDefects(data_carried, [], [], preprocessing_params.defect_mask);
        %else
            [data_masked, preprocessing_params.defect_mask2, defect_centers2, sigmas2] = gaussianMaskDefects(data_carried, preprocessing_params.defect_slice, preprocessing_params.num_defect_type);
        %end
    case 'threshold'
        % Apply threshold-based defect masking
        [data_masked, defect_mask] = thresholdDefects(data_carried, preprocessing_params.defect_slice);
    otherwise
        error('Unknown defect masking method. Choose "gw", "disk", or "threshold".');
end
data_carried = data_masked;

%% 2.6a: Correct streak
% Normalize background 
[data_carried] = normalizeBackgroundToZeroMean3D(data_carried, rangetype, preprocessing_params.slice_normalize);

[data_streakremoved, QPI_nostreaks] = RemoveStreaks(data_carried, 'Direction', 'vertical');
data_carried = data_streakremoved;

%% 2.6b: heal
% Normalize background 
[data_carried] = normalizeBackgroundToZeroMean3D(data_carried, rangetype, preprocessing_params.slice_normalize);

preprocessing_params.heal_direction = input('Enter direction to heal (horizontal/vertical/none): ', 's');

data_streakremoved_healed = heal_streaks(data_carried, preprocessing_params.heal_direction);

%% 2.6c: directional_plane (optional, zero slope at one direction)
% Normalize background 
[data_carried] = normalizeBackgroundToZeroMean3D(data_carried, rangetype, preprocessing_params.slice_normalize);

preprocessing_params.real_space_direction = 'horizontal';
[data_plane, mask] = d3plane_directional(data_carried, preprocessing_params.real_space_direction, 'LineWidth', 5);


%% 2.end: Normalize background 
[Y] = normalizeBackgroundToZeroMean3D(data_carried, rangetype, preprocessing_params.slice_normalize);

%% 3. ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Save the preprocessed data~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
save('ZrSiTe0303_FULL.mat', 'data_original', 'Y', 'data_cropped','data_masked',"energy_range", 'preprocessing_params')

