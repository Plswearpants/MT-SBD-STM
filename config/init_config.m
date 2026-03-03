function cfg = init_config(cfg)
%INIT_CONFIG Initialize/upgrade config schema (no hard-coded defaults).
%   cfg = init_config() creates an empty config struct with the expected
%   nested fields present. This is intended to be called after a checkpoint
%   node is selected/loaded, to ensure any new config fields exist.
%
%   Design rule: this function should NOT hard-code numerical defaults that
%   affect scientific results. Defaults belong in the trunk script PRESETS
%   sections (or are carried from the parent checkpoint).

    if nargin < 1 || isempty(cfg)
        cfg = struct();
    end

    % Top-level groups expected by wrappers/scripts
    cfg = ensureStruct(cfg, 'load');
    cfg = ensureStruct(cfg, 'reference');
    cfg = ensureStruct(cfg, 'sliceRun');
    cfg = ensureStruct(cfg, 'sliceRunPadded');
    cfg = ensureStruct(cfg, 'isolation');
    cfg = ensureStruct(cfg, 'blockInit');
    cfg = ensureStruct(cfg, 'blockRun');
    cfg = ensureStruct(cfg, 'io');

    % Leaf fields (schema only; initialize to [] if missing)
    cfg.load = ensureField(cfg.load, 'data_file', []);
    cfg.load = ensureField(cfg.load, 'smoothing_sigma', []);

    cfg.reference = ensureField(cfg.reference, 'default_ref_slice', []);
    cfg.reference = ensureField(cfg.reference, 'default_num_kernels', []);
    cfg.reference = ensureField(cfg.reference, 'same_size', []);
    cfg.reference = ensureField(cfg.reference, 'kerneltype', []);
    cfg.reference = ensureField(cfg.reference, 'window_type', []);
    cfg.reference = ensureField(cfg.reference, 'square_size', []);

    cfg.sliceRun = ensureField(cfg.sliceRun, 'miniloop_iteration', []);
    cfg.sliceRun = ensureField(cfg.sliceRun, 'outerloop_maxIT', []);
    cfg.sliceRun = ensureField(cfg.sliceRun, 'lambda1', []);
    cfg.sliceRun = ensureField(cfg.sliceRun, 'phase2', []);
    cfg.sliceRun = ensureField(cfg.sliceRun, 'kplus_factor', []);
    cfg.sliceRun = ensureField(cfg.sliceRun, 'lambda2', []);
    cfg.sliceRun = ensureField(cfg.sliceRun, 'nrefine', []);
    cfg.sliceRun = ensureField(cfg.sliceRun, 'signflip', []);
    cfg.sliceRun = ensureField(cfg.sliceRun, 'xpos', []);
    cfg.sliceRun = ensureField(cfg.sliceRun, 'getbias', []);
    cfg.sliceRun = ensureField(cfg.sliceRun, 'Xsolve', []);

    cfg.sliceRunPadded = ensureField(cfg.sliceRunPadded, 'target_size', []);
    cfg.sliceRunPadded = ensureField(cfg.sliceRunPadded, 'miniloop_iteration', []);
    cfg.sliceRunPadded = ensureField(cfg.sliceRunPadded, 'outerloop_maxIT', []);
    cfg.sliceRunPadded = ensureField(cfg.sliceRunPadded, 'lambda1', []);
    cfg.sliceRunPadded = ensureField(cfg.sliceRunPadded, 'phase2', []);
    cfg.sliceRunPadded = ensureField(cfg.sliceRunPadded, 'kplus_factor', []);
    cfg.sliceRunPadded = ensureField(cfg.sliceRunPadded, 'lambda2', []);
    cfg.sliceRunPadded = ensureField(cfg.sliceRunPadded, 'nrefine', []);
    cfg.sliceRunPadded = ensureField(cfg.sliceRunPadded, 'signflip', []);
    cfg.sliceRunPadded = ensureField(cfg.sliceRunPadded, 'xpos', []);
    cfg.sliceRunPadded = ensureField(cfg.sliceRunPadded, 'getbias', []);
    cfg.sliceRunPadded = ensureField(cfg.sliceRunPadded, 'Xsolve', []);

    cfg.isolation = ensureField(cfg.isolation, 'target_kernel_size_type', []);
    cfg.isolation = ensureField(cfg.isolation, 'activation_threshold_divisor', []);

    cfg.blockInit = ensureField(cfg.blockInit, 'use_matrix', []);
    cfg.blockInit = ensureField(cfg.blockInit, 'change_size', []);

    cfg.blockRun = ensureField(cfg.blockRun, 'trusted_ratio_threshold_default', []);
    cfg.blockRun = ensureField(cfg.blockRun, 'use_default_manual_trusted_slices', []);
    cfg.blockRun = ensureField(cfg.blockRun, 'show_trusted_plot', []);
    cfg.blockRun = ensureField(cfg.blockRun, 'miniloop_iteration', []);
    cfg.blockRun = ensureField(cfg.blockRun, 'outerloop_maxIT', []);
    cfg.blockRun = ensureField(cfg.blockRun, 'lambda1_base', []);
    cfg.blockRun = ensureField(cfg.blockRun, 'phase2', []);
    cfg.blockRun = ensureField(cfg.blockRun, 'kplus_factor', []);
    cfg.blockRun = ensureField(cfg.blockRun, 'lambda2', []);
    cfg.blockRun = ensureField(cfg.blockRun, 'nrefine', []);
    cfg.blockRun = ensureField(cfg.blockRun, 'signflip', []);
    cfg.blockRun = ensureField(cfg.blockRun, 'xpos', []);
    cfg.blockRun = ensureField(cfg.blockRun, 'getbias', []);
    cfg.blockRun = ensureField(cfg.blockRun, 'Xsolve', []);
    cfg.blockRun = ensureField(cfg.blockRun, 'use_Xregulated', []);
    cfg.blockRun = ensureField(cfg.blockRun, 'allow_custom_update_order', []);

    cfg.io = ensureField(cfg.io, 'allslice_output_file', []);
end

function s = ensureStruct(s, fieldName)
    if ~isfield(s, fieldName) || ~isstruct(s.(fieldName))
        s.(fieldName) = struct();
    end
end

function s = ensureField(s, fieldName, defaultValue)
    if ~isfield(s, fieldName)
        s.(fieldName) = defaultValue;
    end
end

