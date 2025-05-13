QPI_original = qpiCalculate(data_streakremoved_healed);
[masks, masks_combined, comment] = maskDirectional(QPI_original);

%% 
list = 1:size(masks,3);
slist= list(1:20:end);

p = 100;
masked_data = zeros(size(masks,1),size(masks,2),size(slist,2)); 
%%
for s = list(1:20:end)
    masked_data(:,:,find(slist==s))= QPI_original(:,:,p).* masks(:,:,s);
end

masked_data_trim = trim_zeros(masked_data);
%% plot 
X = 1:size(masked_data_trim,1); 
height = mean(masked_data_trim,'all');

% Create a colormap with 'num_lines' colors
cmap = parula(size(masked_data_trim,3));  % can also use 'jet', 'turbo', 'gray', etc.

figure; hold on;
for i = 1:size(masked_data_trim,3)
    height = height * 1.1;
    y = height + reshape(masked_data_trim(:,:,i), [], 1);  % forces column vector
    plot(X,y,'Color', cmap(i,:), 'LineWidth', 2)
end
hold off;

colorbar;  % optional, shows color scale



%%
function A_trimmed = trim_zeros(A)
    [nx, ny, nz] = size(A);

    % First pass: find bounding boxes for each slice
    bounds = zeros(nz, 4); % [row_min, row_max, col_min, col_max]
    for k = 1:nz
        [rows, cols] = find(A(:, :, k));
        if ~isempty(rows)
            bounds(k, :) = [min(rows), max(rows), min(cols), max(cols)];
        else
            bounds(k, :) = [1, 1, 1, 1];  % dummy 1x1 zero patch if slice is empty
        end
    end

    % Determine max output size
    row_sizes = bounds(:,2) - bounds(:,1) + 1;
    col_sizes = bounds(:,4) - bounds(:,3) + 1;
    max_row = max(row_sizes);
    max_col = max(col_sizes);

    % Preallocate output
    A_trimmed = zeros(max_row, max_col, nz);

    % Second pass: extract and pad each slice
    for k = 1:nz
        r1 = bounds(k,1); r2 = bounds(k,2);
        c1 = bounds(k,3); c2 = bounds(k,4);
        slice = A(r1:r2, c1:c2, k);

        [h, w] = size(slice);
        A_trimmed(1:h, 1:w, k) = slice;
    end
end

