function combined_metrics = combine_metrics_for_plotting(dataset_metrics_array)
% Combines multiple dataset_metrics structures into one, averaging overlapping points
% Output is compatible with plot_defects_snr_kernel_similarity

% Collect all unique parameter values
all_snr = [];
all_theta = [];
all_nobs = [];
for i = 1:length(dataset_metrics_array)
    m = dataset_metrics_array{i};
    all_snr = [all_snr; m.SNR_values(:)];
    all_theta = [all_theta; m.theta_cap_values(:)];
    all_nobs = [all_nobs; m.Nobs_values(:)];
end
unique_snr = unique(all_snr);
unique_theta = unique(all_theta);
unique_nobs = unique(all_nobs);

% Create combined arrays
sz = [length(unique_snr), length(unique_theta), length(unique_nobs)];
combined_kernel_quality = nan(sz);
combined_X0 = cell(sz);
point_counts = zeros(sz);

% For each run, map its metrics into the combined arrays
for run_idx = 1:length(dataset_metrics_array)
    m = dataset_metrics_array{run_idx};
    for i = 1:length(m.SNR_values)
        for j = 1:length(m.theta_cap_values)
            for k = 1:length(m.Nobs_values)
                % Find indices in combined arrays
                [~, snr_idx] = min(abs(unique_snr - m.SNR_values(i)));
                [~, theta_idx] = min(abs(unique_theta - m.theta_cap_values(j)));
                [~, nobs_idx] = min(abs(unique_nobs - m.Nobs_values(k)));
                % Get kernel quality
                kq = m.kernel_quality_final(i,j,k);
                X0 = m.X0{i,j,k};
                if isnan(kq) || isempty(X0)
                    continue;
                end
                % Average overlapping points
                if isnan(combined_kernel_quality(snr_idx, theta_idx, nobs_idx))
                    combined_kernel_quality(snr_idx, theta_idx, nobs_idx) = kq;
                    combined_X0{snr_idx, theta_idx, nobs_idx} = X0;
                    point_counts(snr_idx, theta_idx, nobs_idx) = 1;
                else
                    % Average kernel quality
                    total = combined_kernel_quality(snr_idx, theta_idx, nobs_idx) * point_counts(snr_idx, theta_idx, nobs_idx);
                    total = total + kq;
                    point_counts(snr_idx, theta_idx, nobs_idx) = point_counts(snr_idx, theta_idx, nobs_idx) + 1;
                    combined_kernel_quality(snr_idx, theta_idx, nobs_idx) = total / point_counts(snr_idx, theta_idx, nobs_idx);
                    % For X0, just keep the first non-empty one (or you could average them if needed)
                end
            end
        end
    end
end

combined_metrics = struct();
combined_metrics.kernel_quality_final = combined_kernel_quality;
combined_metrics.X0 = combined_X0;
combined_metrics.SNR_values = unique_snr;
end 