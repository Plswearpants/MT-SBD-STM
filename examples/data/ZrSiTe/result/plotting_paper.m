figure; d3gridDisplay(data.real.Y, 'dynamic');

%% Figure 1 
% def params
obs = Y_used;
chosen_slice = 80; 
num_kernels = 2;

%% produce cropping 
square_size = [80,80];
kernel_sizes = repmat(square_size,[num_kernels,1]);
[A1_ref, ~] = initialize_kernels(obs(:,:,chosen_slice), num_kernels, kernel_sizes, 'selected', '');

%% Make plot
fig1_suba = obs(:,:,chosen_slice);

figure;
title('subplot a')
imagesc(fig1_suba); axis square; colormap('gray');
set(gca,'XTick', [], 'YTick', [])

fig1_subb = qpiCalculate(obs(:,:,chosen_slice));
figure; 
title('subplot b')
imagesc(fig1_subb); axis square; colormap('invgray');
clim = clipEdgeIntensity(fig1_subb, 1);
set(gca,'XTick', [], 'YTick', [], 'CLim', clim)

figure; 
title('cropping')
ax1 = subplot(2,2,1);
imagesc(A1_ref{1});axis square; colormap(ax1,'gray');
set(ax1, 'XTick', [], 'YTick', []);
ax2 = subplot(2,2,2);
imagesc(qpiCalculate(A1_ref{1}));axis square; colormap(ax2,'invgray');
set(ax2, 'XTick', [], 'YTick', [], 'Clim', clipEdgeIntensity(qpiCalculate(A1_ref{1}), 1))
ax3 = subplot(2,2,3);
imagesc(A1_ref{2});axis square; colormap(ax3, 'gray');
set(ax3, 'XTick', [], 'YTick', [])
ax4 = subplot(2,2,4);
imagesc(qpiCalculate(A1_ref{2}));axis square; colormap(ax4,'invgray');
set(ax4, 'XTick', [], 'YTick', [], 'Clim', clipEdgeIntensity(qpiCalculate(A1_ref{2}), 1))

%% Figure 5 
angle = 122.5741;
Xout_A1 = Xout_A1(:,:,1);
Xout_A2 = Xout_A2(:,:,1);
target_size =[242,242];
radius = 242;
%% Create padded and symmetrized QPI_out
% pad QPI
pad_QPI_each = zeros([num_kernels,size(Y_used)]);
for i = 1:num_kernels
    pad_QPI_each(i,:,:,:) = qpiCalculate(Aout_Full_energy{i,1},size(Y_used,[1,2]));
end

% Symmetrize the data
[Y_rec_symm, Y_rec_symm_45,tform0, angle]   = Symmetrizing2(squeeze(qpiCalculate(Y_used)),'',angle, 'default');
[A1_rec_symm, A1_rec_symm_45,tform1, angle] = Symmetrizing2(squeeze(pad_QPI_each(1,:,:,:)),'',angle, 'default');
[A2_rec_symm, A2_rec_symm_45,tform2, angle] = Symmetrizing2(squeeze(pad_QPI_each(2,:,:,:)),'',angle, 'default');

% Crop symmetrized output to match QPI_rec x-y size
%target_size = size(Y_used);
Y_cropped = centerCropToTargetSize(Y_rec_symm_45, target_size);
A1_cropped = centerCropToTargetSize(A1_rec_symm_45, target_size);
A2_cropped = centerCropToTargetSize(A2_rec_symm_45, target_size);

%% Compare A1,A2 and Y in q-space
QPI_Y_norm = normalizeForDisplay(Y_cropped);
A1_cropped_norm = normalizeForDisplay(A1_cropped);
A2_cropped_norm = normalizeForDisplay(A2_cropped);
comparison = [QPI_Y_norm, A1_cropped_norm, A2_cropped_norm];
figure; d3gridDisplay(comparison, 'dynamic', -1)

% [X1_sum_points, X1_num_points] = collectAboveThreshold(Xout_A1, 0.1);
% [X2_sum_points, X2_num_points] = collectAboveThreshold(Xout_A2, 0.1);
% comparison_no_norm = [Y_cropped, X1_sum_points*A1_cropped, X2_sum_points*A2_cropped];
% figure; d3gridDisplay(comparison_no_norm, 'dynamic', -1)

%% Plot rotating slices 

% Create rotational slices
[A1_rotational_slices, A1_slice_angles, ~] = rotationalslices(A1_cropped, 'global', 1, 121);
[A2_rotational_slices, A2_slice_angles, ~] = rotationalslices(A2_cropped, 'global', 1, 121);
%% -------------------------------------------------------------------------
% Helper: clip intensity to (n, 100-n) percentiles for better readability
% Usage: clim = clipEdgeIntensity(Z, n);  clim(ax, clim);  % or caxis(clim)
function clim = clipEdgeIntensity(Z, n)
    if nargin < 2, n = 2; end
    n = max(0, min(50, n));
    clim = [prctile(Z(:), n), prctile(Z(:), 100 - n)];
end

%% -------------------------------------------------------------------------
% Helper: set second axes intensity window to match the first (reference)
% Usage: unifyIntensityWindow(ax_ref, ax_to_update)
function unifyIntensityWindow(ax_ref, ax_to_update)
    set(ax_to_update, 'CLim', get(ax_ref, 'CLim'));
end

%% -------------------------------------------------------------------------
% Helper: crop center region of 3D data to target [rows, cols]
% Usage: out = centerCropToTargetSize(data3d, [m, n])
function data_cropped = centerCropToTargetSize(data_tobe_cropped, target_size)
    validateattributes(data_tobe_cropped, {'numeric'}, {'nonempty'});
    validateattributes(target_size, {'numeric'}, {'vector','numel',2,'integer','positive','finite','real'});

    size_data = [size(data_tobe_cropped,1), size(data_tobe_cropped,2)];
    if any(target_size > size_data)
        error('centerCropToTargetSize:TargetTooLarge', ...
            'target_size [%d %d] exceeds data size [%d %d].', ...
            target_size(1), target_size(2), size_data(1), size_data(2));
    end

    x_range = [1 + floor((size_data(1)-target_size(1))/2), floor((size_data(1)+target_size(1))/2)];
    y_range = [1 + floor((size_data(2)-target_size(2))/2), floor((size_data(2)+target_size(2))/2)];
    data_cropped = data_tobe_cropped(x_range(1):x_range(2), y_range(1):y_range(2), :);
end

%% -------------------------------------------------------------------------
% Helper: normalize numeric array to [0, 1] for visualization
% Usage: out = normalizeForDisplay(data)
function data_norm = normalizeForDisplay(data_in)
    validateattributes(data_in, {'numeric'}, {'nonempty', 'real', 'finite'});
    min_val = min(data_in(:));
    max_val = max(data_in(:));
    range_val = max_val - min_val;
    if range_val == 0
        data_norm = zeros(size(data_in), 'like', data_in);
    else
        data_norm = (data_in - min_val) ./ range_val;
    end
end

%% -------------------------------------------------------------------------
% Helper: sum and count points in 2D map above threshold
% Usage: [sum_points, num_points] = collectAboveThreshold(Xout, 0.1)
function [sum_points, num_points] = collectAboveThreshold(Xout, threshold)
    if nargin < 2
        threshold = 0.1;
    end
    validateattributes(Xout, {'numeric'}, {'2d', 'nonempty', 'real'});
    validateattributes(threshold, {'numeric'}, {'scalar', 'real', 'finite'});

    selected_points = Xout(Xout > threshold);
    sum_points = sum(selected_points(:));
    num_points = numel(selected_points);
end