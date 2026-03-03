function [log, data, params, meta, cfg] = loadRealDataset(log, data, params, meta, cfg)
%LOADREALDATASET Load real STM/QPI .3ds data into structured form.
%
%   [log, data, params, meta, cfg] = loadRealDataset(log, data, params, meta, cfg)
%
%   This is a thin wrapper around the legacy Block 1 in
%   scripts/MTSBD_block_realdata1.m. It:
%       - Ensures cfg.load.data_file and cfg.load.smoothing_sigma are set
%       - Calls load3dsall to load the .3ds file
%       - Builds basic energy axis and dimensions
%       - Stores results under data.real and params.real
%
%   Project/checkpoint creation and logging are handled at the script level
%   (see scripts/run_real_data.m and docs/script_standardization.md).
%
%   Inputs:
%       log, data, params, meta, cfg  - pipeline structs (may be empty)
%
%   Outputs:
%       Updated log, data.real, params.real, meta, cfg

    arguments
        log  struct
        data struct
        params struct
        meta struct
        cfg  struct
    end

    % Ensure load sub-struct exists
    if ~isfield(cfg, "load")
        error('cfg.load must be defined (see get_config).');
    end

    % If data_file is missing or empty, allow interactive selection
    if ~isfield(cfg.load, "data_file") || isempty(cfg.load.data_file)
        [f, p] = uigetfile({'*.3ds','Nanonis 3DS files (*.3ds)'}, ...
                           'Select real-data .3ds file');
        if isequal(f, 0)
            error('No .3ds file selected.');
        end
        cfg.load.data_file = fullfile(p, f);
    end

    if ~isfield(cfg.load, "smoothing_sigma") || isempty(cfg.load.smoothing_sigma)
        cfg.load.smoothing_sigma = 10;
    end

    % Call legacy loader (Block 1 logic)
    [header, par, I, dIdV, LockindIdV, bias, midV, QPI, LockinQPI] = ...
        load3dsall(cfg.load.data_file, cfg.load.smoothing_sigma);

    xsize  = header.grid_dim(1);
    ysize  = header.grid_dim(2);
    elayer = header.points;
    estart = par(1);
    eend   = par(2);
    energy_range = linspace(estart, eend, elayer);

    data_original = dIdV;
    num_slices    = size(data_original, 3);
    spatial       = size(data_original, 1);

    % Populate data.real
    if ~isfield(data, "real"); data.real = struct(); end

    data.real.header        = header;
    data.real.par           = par;
    data.real.I             = I;
    data.real.dIdV          = dIdV;
    data.real.LockindIdV    = LockindIdV;
    data.real.bias          = bias;
    data.real.midV          = midV;
    data.real.QPI           = QPI;
    data.real.LockinQPI     = LockinQPI;
    data.real.data_original = data_original;
    data.real.energy_range  = energy_range;
    data.real.xsize         = xsize;
    data.real.ysize         = ysize;
    data.real.elayer        = elayer;

    % Basic params.real fields (dimensions and energy)
    if ~isfield(params, "real"); params.real = struct(); end

    params.real.num_slices   = num_slices;
    params.real.spatial_size = spatial;
    params.real.energy_start = estart;
    params.real.energy_end   = eend;

    % Minimal metadata for future checkpoint logic
    if ~isfield(meta, "raw_path_project") || isempty(meta.raw_path_project)
        % For now, just record the source file; project-local copying will
        % be handled in the checkpoint/IO layer.
        meta.raw_path_project = cfg.load.data_file;
    end

end

