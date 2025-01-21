function IDX = calculate_overlap_regions(Xiter, k)
    % this function calculate the overlap regions between the activation regions
    % Xiter: the activation regions
    % k: the size of the kernel
    % IDX: the overlap region, size = [h, w, num_kernels, num_kernels]
    %      IDX(:,:,i,j) shows where kernel i overlaps with kernel j

    [h, w, num_kernels] = size(Xiter);
    IDX = false(h, w, num_kernels, num_kernels);
    Mask = false(size(Xiter));
    
    % create square mask around every activation region
    for n = 1:num_kernels
        Mask(:,:,n) = conv2(Xiter(:,:,n), ones(k(n,1), k(n,2)), 'same') > 0;
    end
    
    % calculate individual overlaps between each pair of kernels
    for n = 1:num_kernels
        for m = 1:num_kernels
            if n ~= m
                IDX(:,:,n,m) = Mask(:,:,n) & Mask(:,:,m);
            end
        end
    end
end
