%% Post processing:
% Aout_All is a n by 1 cell, each cell is a 3D matrix. Now we start post processing output kernels.
% first we mask out the kernel center with truncated gaussian window, following similar procedure as in Block 2.5.
% specify which kernel type to process
kernel_type_to_process = [3,4];
for k = kernel_type_to_process
    f1=figure;
    d3gridDisplay(Aout_ALL{k},'dynamic');
    preprocessing_params.defect_slice = input('Enter defect slice number: ');
    preprocessing_params.num_defect_type = input('enter how many types of defects to mask: ');
    close(f1);
    [Aout_ALL{k}, preprocessing_params.defect_mask, defect_centers, sigmas] = gaussianMaskDefects(Aout_ALL{k}, preprocessing_params.defect_slice, preprocessing_params.num_defect_type);
end

%% normalize method 1: each slice for individual type of defect
Aout_show_Full = [];
qpi_show_Full = [];
for i = 1: size(Aout_ALL,1)
    pp = Aout_ALL{i};
    qq = Aout_ALL{i};
    for j = 1: size(pp,3)
        pp(:,:,j) = mat2gray(pp(:,:,j));
        qq(:,:,j) = 1-mat2gray(qpiCalculate(qq(:,:,j)),[0,1]);
    end
    Aout_show_Full = [Aout_show_Full,pp];
    qpi_show_Full = [qpi_show_Full,qq];
end
ALL_show_norm = [Aout_show_Full;qpi_show_Full];

%% normalize method 2: each slice for all kernels
Aout_show_Full = [];
qpi_show_Full = [];
for i = 1: size(Aout_ALL,1)
    Aout_show_Full = [Aout_show_Full,Aout_ALL{i}];
    qpi_show_Full = [qpi_show_Full,qpiCalculate(Aout_ALL{i})];
end

for i = 1: size(Aout_show_Full,3)
    Aout_show_Full(:,:,i) = mat2gray(Aout_show_Full(:,:,i));
    qpi_show_Full(:,:,i) = 1-mat2gray(qpi_show_Full(:,:,i),[0,1]);
end
ALL_show_norm = [Aout_show_Full;qpi_show_Full];
%%
figure; 
d3gridDisplay(ALL_show_norm, 'dynamic');

%% Show FULL
Aout_show_Full = [];
for i = 1: size(Aout_ALL,1)
    Aout_show_Full = [Aout_show_Full,Aout_ALL{i}];
end

figure;
d3gridDisplay(Aout_show_Full, 'dynamic')

%%
qpi_show_Full = [];
for i = 1: size(Aout_ALL,1)
    qpi_show_Full = [qpi_show_Full,qpiCalculate(Aout_ALL{i})];
end
figure;
d3gridDisplay(qpi_show_Full, 'dynamic',-1)

%%
Aout_show_norm = Aout_show_Full;
qpi_show_norm = qpi_show_Full;
for i = 1: 200
    Aout_show_norm(:,:,i) = mat2gray(Aout_show_Full(:,:,i));
    qpi_show_norm(:,:,i) = abs(mat2gray(qpi_show_Full(:,:,i),[0,1])-1);
end
ALL_show_norm = [Aout_show_norm;qpi_show_norm];
figure; 
d3gridDisplay(ALL_show_norm, 'dynamic');

%% write the video 
gridVideoWriter(rot90(ALL_show_norm), V, 'dynamic', 50, 'gray', 0, [800 800]);


%% convert Aout_ALL to cell format
%[num_slices, num_kernels] = size(bout_ALL);
Aout_ALL_cell = cell(num_slices, num_kernels);
for s = 1:num_slices
    for k = 1:num_kernels
        Aout_ALL_cell{s,k} = Aout_Full_energy{k}(:,:,s);
    end
end

%% Create reconstruction for each kernel type
Y_rec_each = zeros([num_kernels,size(Y_used)]);
for i = 1:size(Y_used,3)
    for k = 1:num_kernels
        %Y_rec_each(k,:,:,i) = convfft2(Aout_ALL_cell{i,k}, Xout_ALL(:,:,k)) + bout_ALL(i,k);
        Y_rec_each(k,:,:,i) = convfft2(Aout_ALL_cell{i,k}, Xout_ALL(:,:,k));
    end
end
% create fft of Y_rec_each
FT_QPI_Y_rec_each = zeros([num_kernels,size(Y_used)]);
for k = 1:num_kernels
    FT_QPI_Y_rec_each(k,:,:,:) = qpiCalculate(squeeze(Y_rec_each(k,:,:,:)));
end

%% Normalize and combine Y_rec_each and its FT-QPI using method 2
% Reshape Y_rec_each and FT_QPI_Y_rec_each to combine all kernels
Y_rec_show_Full = [];
qpi_Y_rec_show_Full = [];
for k = 1:num_kernels
    Y_rec_show_Full = [Y_rec_show_Full, squeeze(Y_rec_each(k,:,:,:))];
    qpi_Y_rec_show_Full = [qpi_Y_rec_show_Full, squeeze(FT_QPI_Y_rec_each(k,:,:,:))];
end

% Normalize each slice across all kernels
for i = 1:size(Y_rec_show_Full,3)
    Y_rec_show_Full(:,:,i) = mat2gray(Y_rec_show_Full(:,:,i));
    qpi_Y_rec_show_Full(:,:,i) = 1-mat2gray(qpi_Y_rec_show_Full(:,:,i),[0,1]);
end

% Combine normalized reconstructions and their FT-QPI
Y_rec_ALL_show_norm = [Y_rec_show_Full; qpi_Y_rec_show_Full];

%% Display the normalized and combined results
figure;
d3gridDisplay(Y_rec_ALL_show_norm, 'dynamic');
title('Normalized Y_{rec}_-{each} and FT-QPI combined');