qpi_Aout_pad = zeros(477,5*477,51);

for i = 1:5
    qpi_Aout_pad(:,1+(i-1)*477:i*477,:) = qpiCalculate(Aout_ALL{i},477);
end

% Normalize each slice across all kernels
for i = 1:size(qpi_Aout_pad,3)
    qpi_Aout_pad(:,:,i) = 1-mat2gray(qpi_Aout_pad(:,:,i),[0,1]);
end

Y_rec_ALL_show_norm= [Y_rec_ALL_show_norm;qpi_Aout_pad];

%%

qpi_Aout_resize = zeros(477,5*477,51);

for i = 1:5
    qpi_Aout_resize(:,1+(i-1)*477:i*477,:) = imresize(qpiCalculate(Aout_ALL{i},80),[477,477]);
end

% Normalize each slice across all kernels
for i = 1:size(qpi_Aout_resize,3)
    qpi_Aout_resize(:,:,i) = 1-mat2gray(qpi_Aout_resize(:,:,i),[0,1]);
end

Y_rec_ALL_show_norm= [Y_rec_ALL_show_norm;qpi_Aout_resize];

%%
padding_l = 250;
qpi_Aout_pad_any = zeros(padding_l,5*padding_l,51);

for i = 1:5
    qpi_Aout_pad_any(:,1+(i-1)*padding_l:i*padding_l,:) = qpiCalculate(Aout_ALL{i},padding_l);
end

% Normalize each slice across all kernels
for i = 1:size(qpi_Aout_pad_any,3)
    qpi_Aout_pad_any(:,:,i) = 1-mat2gray(qpi_Aout_pad_any(:,:,i),[0,1]);
end

%%
center = [1241,249];
side = 120;
nn = [combined_data(249-120:249+120,1241-120:1241+120,:);combined_data(745-120:745+120,1241-120:1241+120,:)];