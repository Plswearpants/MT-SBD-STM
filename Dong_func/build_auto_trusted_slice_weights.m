function [slice_weights, details] = build_auto_trusted_slice_weights(A1_all_matrix, noise_var, threshold, show_plot, manual_trusted_slices)
%BUILD_AUTO_TRUSTED_SLICE_WEIGHTS Build binary trusted-slice weights per kernel.
%   Inputs:
%       A1_all_matrix : cell array, each cell is [h x w x num_slices] kernel stack
%       noise_var     : scalar or [num_slices x 1] noise variance
%       threshold     : trusted ratio threshold
%       show_plot     : whether to show ratio-vs-slice visualization
%       manual_trusted_slices : cell array of manual trusted slice indices (optional)
%   Rule:
%       ratio(i,j) = std(A1_all_matrix{j}(:,:,i)) / sqrt(noise_var(i))
%       trusted if ratio(i,j) > threshold

    if nargin < 4 || isempty(show_plot)
        show_plot = true;
    end
    if nargin < 5 || isempty(manual_trusted_slices)
        manual_trusted_slices = {};
    end

    num_kernels = numel(A1_all_matrix);
    num_slices = size(A1_all_matrix{1}, 3);
    eps0 = 1e-12;

    if isscalar(noise_var)
        noise_var = repmat(noise_var, [num_slices, 1]);
    else
        noise_var = noise_var(:);
    end
    if numel(noise_var) ~= num_slices
        error('noise_var must be scalar or one value per slice.');
    end

    noise_std = sqrt(max(double(noise_var), eps0));
    ratio = zeros(num_slices, num_kernels);
    slice_weights = zeros(num_slices, num_kernels);
    trusted_slices = cell(1, num_kernels);
    trusted_counts = zeros(1, num_kernels);

    for k = 1:num_kernels
        Ak = A1_all_matrix{k};
        if size(Ak, 3) ~= num_slices
            error('All kernel stacks in A1_all_matrix must have the same num_slices.');
        end

        for s = 1:num_slices
            curr_slice = double(Ak(:,:,s));
            ratio(s,k) = std(curr_slice(:), 0) / noise_std(s);
        end

        idx = find(ratio(:,k) > threshold);
        if isempty(idx)
            % Ensure at least one trusted slice: keep the max-ratio slice.
            [~, idx_max] = max(ratio(:,k));
            idx = idx_max;
        end

        slice_weights(idx, k) = 1;
        trusted_slices{k} = idx(:).';
        trusted_counts(k) = numel(idx);
    end

    if show_plot
        figure('Name', 'Trusted Slice Ratio by Kernel');
        hold on;
        h = gobjects(num_kernels,1);
        for k = 1:num_kernels
            h(k) = scatter(1:num_slices, ratio(:,k), 18, 'filled');
            if numel(manual_trusted_slices) >= k && ~isempty(manual_trusted_slices{k})
                idxm = unique(round(manual_trusted_slices{k}(:).'));
                idxm = idxm(idxm >= 1 & idxm <= num_slices);
                if ~isempty(idxm)
                    c = h(k).CData;
                    scatter(idxm, ratio(idxm,k), 60, c, 'o', 'LineWidth', 1.2);
                end
            end
        end
        yline(threshold, 'k--', 'LineWidth', 1.5);
        xlabel('Slice index');
        ylabel('std(kernel slice) / std(noise)');
        title('Trusted-slice std ratio vs slice index');
        legend_labels = arrayfun(@(k) sprintf('Kernel %d', k), 1:num_kernels, 'UniformOutput', false);
        h_manual = scatter(nan, nan, 60, 'o', 'k', 'LineWidth', 1.2);
        h_thr = plot(nan, nan, 'k--', 'LineWidth', 1.5);
        legend([h; h_manual; h_thr], [legend_labels, {'Manual trusted slices', 'Threshold'}], 'Location', 'best');
        grid on;
        hold off;
    end

    details = struct();
    details.threshold = threshold;
    details.ratio = ratio;
    details.trusted_slices = trusted_slices;
    details.trusted_counts = trusted_counts;
end