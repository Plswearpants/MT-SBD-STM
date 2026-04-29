function apply_defect_density_tick_style(ax, theta_values)
% Standardize defect-density ticks to clean log-decade labels.
% Example labels: 10^{-3}, 10^{-2}, 10^{-1}
    if nargin < 1 || isempty(ax) || ~isgraphics(ax)
        ax = gca;
    end
    if nargin < 2 || isempty(theta_values)
        return;
    end

    theta_values = theta_values(:)';
    theta_values = theta_values(isfinite(theta_values) & theta_values > 0);
    if isempty(theta_values)
        return;
    end
    theta_values = unique(theta_values);

    tmin = min(theta_values);
    tmax = max(theta_values);
    pmin = floor(log10(tmin));
    pmax = ceil(log10(tmax));
    major_ticks = 10.^(pmin:pmax);

    in_range = major_ticks >= tmin * (1 - 1e-12) & major_ticks <= tmax * (1 + 1e-12);
    major_ticks = major_ticks(in_range);
    if isempty(major_ticks)
        major_ticks = theta_values([1 end]);
    end

    % Explicitly define minor ticks (2..9 per decade) so the spacing is
    % consistent and visible between adjacent powers of ten.
    minor_ticks = [];
    if pmax >= pmin
        decade_grid = 10.^(pmin:pmax);
        for d = 1:numel(decade_grid)
            minor_ticks = [minor_ticks, (2:9) * decade_grid(d)]; %#ok<AGROW>
        end
        in_minor_range = minor_ticks >= tmin * (1 - 1e-12) & minor_ticks <= tmax * (1 + 1e-12);
        minor_ticks = unique(minor_ticks(in_minor_range));
    end

    ax.XScale = 'log';
    ax.XTick = major_ticks;
    ax.XTickLabel = arrayfun(@(v) sprintf('10^{%d}', round(log10(v))), major_ticks, 'UniformOutput', false);
    ax.XTickLabelRotation = 0;
    ax.XMinorTick = 'on';
    if isprop(ax, 'XAxis') && isprop(ax.XAxis, 'MinorTickValues')
        try
            ax.XAxis.MinorTickValues = minor_ticks;
        catch
            % Keep default automatic minor ticks if explicit assignment fails.
        end
    end
end
