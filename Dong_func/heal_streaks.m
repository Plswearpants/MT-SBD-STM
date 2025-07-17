function data_streakremoved_healed = heal_streaks(data_carried, direction)
%HEAL_STREAKS Heal the midlines of data_streakremoved in the specified direction.
%   data_streakremoved_healed = heal_streaks(data_carried, direction)
%   direction: 'horizontal', 'vertical', or 'none'

    data_streakremoved_healed = data_carried; % Default if no healing is applied
    switch lower(direction)
        case 'horizontal'
            QPI = zeros(size(data_carried));
            mid_row = round(size(QPI,1)/2)+1;
            for i = 1:size(data_carried,3)
                QPI(:,:,i) = fftshift(fft2(data_carried(:,:,i)));
                % First heal the middle point using its neighbors
                QPI(mid_row, :, i) = mean([QPI(mid_row-1, :, i); QPI(mid_row+1, :, i)], 1);
                % Then heal mid-1 using its neighbors
                QPI(mid_row-1, :, i) = mean([QPI(mid_row-2, :, i); QPI(mid_row, :, i)], 1);
                % Finally heal mid+1 using its neighbors
                QPI(mid_row+1, :, i) = mean([QPI(mid_row, :, i); QPI(mid_row+2, :, i)], 1);
                data_streakremoved_healed(:,:,i) = real(ifft2(ifftshift(QPI(:,:,i))));
            end
        case 'vertical'
            QPI = zeros(size(data_carried));
            mid_col = round(size(QPI,2)/2);
            for i = 1:size(data_carried,3)
                QPI(:,:,i) = fftshift(fft2(data_carried(:,:,i)));
                % First heal the middle point using its neighbors
                QPI(:, mid_col, i) = mean([QPI(:, mid_col-1, i), QPI(:, mid_col+1, i)], 2);
                % Then heal mid-1 using its neighbors
                QPI(:, mid_col-1, i) = mean([QPI(:, mid_col-2, i), QPI(:, mid_col, i)], 2);
                % Finally heal mid+1 using its neighbors
                QPI(:, mid_col+1, i) = mean([QPI(:, mid_col, i), QPI(:, mid_col+2, i)], 2);
                data_streakremoved_healed(:,:,i) = real(ifft2(ifftshift(QPI(:,:,i))));
            end
        case 'none'
            % No healing needed, use data as is
            disp('No healing applied');
        otherwise
            error('Unknown direction for healing: %s', direction);
    end
end 