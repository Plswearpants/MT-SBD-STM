function [Aout, Xout, bout, extras] = MTSBD_all_slice(Y, k, params, dispfun, kernel_initialguess, Max_iteration, maxIT)
%MTSBD_ALL_SLICE Backward-compatible wrapper for retired implementation.
%
% The legacy implementation has been retired in favor of
% `MTSBD_all_slice_modified`. The pre-retirement body is archived at
% `history/MTSBD_all_slice_legacy_pre_retirement.m`.

    warning(['MTSBD_all_slice is retired and now forwards to ', ...
        'MTSBD_all_slice_modified. Update call sites to use ', ...
        'MTSBD_all_slice_modified directly.']);

    [Aout, Xout, bout, extras] = MTSBD_all_slice_modified( ...
        Y, k, params, dispfun, kernel_initialguess, Max_iteration, maxIT);
end