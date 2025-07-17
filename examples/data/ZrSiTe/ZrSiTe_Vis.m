%% Processing script for ZrSiTe dataset(Post run visualization)


%% How well each Xout aligned

% combine Xout for 3 runs for each kernel
Xout_A1 = cat(3, Xout_ALL1(:,:,1), Xout_ALL2(:,:,1), Xout_ALL3(:,:,1));
Xout_A2 = cat(3, Xout_ALL1(:,:,2), Xout_ALL2(:,:,2), Xout_ALL3(:,:,2));

%% Basic Gaussian window
kernel_size = [80,80;80,80;80,80];

sim_matrix_A1 = cosineSimilarityMatrix(Xout_A1, 'kernel_size', kernel_size);
sim_matrix_A2 = cosineSimilarityMatrix(Xout_A2, 'kernel_size', kernel_size);
