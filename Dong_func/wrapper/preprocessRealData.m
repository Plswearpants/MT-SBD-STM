function [log, data, params, meta, cfg] = preprocessRealData(log, data, params, meta, cfg)
%PREPROCESSREALDATA Preprocess real STM/QPI data (Block 2 logic).
%
%   Supports interactive (record choices + masks) and non-interactive replay.
%   When params.preprocessing.interactive is false, loads recipe from
%   data.real.preprocessing and uses steps.flags; no prompts/UI.

    arguments
        log  struct
        data struct
        params struct
        meta struct
        cfg  struct
    end

    if ~isfield(data, "real") || ~isfield(data.real, "data_original")
        error('preprocessRealData: data.real.data_original is missing. Run loadRealDataset first.');
    end
    if ~isfield(data.real, "energy_range")
        error('preprocessRealData: data.real.energy_range is missing.');
    end

    data_original = data.real.data_original;
    energy_range  = data.real.energy_range;

    if ~isfield(params, "preprocessing")
        params.preprocessing = struct();
    end

    step_names = { ...
        'bragg_remove', 'crop', 'slice_select', 'manual_streak', ...
        'auto_streak', 'defect_mask', 'streak_correct', 'heal', 'directional_plane' ...
    };
    nSteps = numel(step_names);

    interactive = true;
    if isfield(params.preprocessing, 'interactive')
        interactive = params.preprocessing.interactive;
    end
    if ~isfield(params.preprocessing, 'save_checkpoint')
        params.preprocessing.save_checkpoint = true;
    end

    % ----- Non-interactive: load recipe and validate -----
    if ~interactive
        if ~isfield(data, 'real') || ~isfield(data.real, 'preprocessing') || isempty(data.real.preprocessing)
            error('preprocessRealData: non-interactive mode requires data.real.preprocessing (load a checkpoint first).');
        end
        preprocessing_params = data.real.preprocessing;
        if ~isfield(preprocessing_params, 'steps') || ~isfield(preprocessing_params.steps, 'flags')
            error('preprocessRealData: data.real.preprocessing.steps.flags is missing.');
        end
        flags = preprocessing_params.steps.flags;
        if numel(flags) ~= nSteps
            error('preprocessRealData: steps.flags length (%d) does not match step_names (%d).', numel(flags), nSteps);
        end
        % Require slice_normalize for final normalization
        if ~isfield(preprocessing_params, 'slice_normalize') || isempty(preprocessing_params.slice_normalize)
            error('preprocessRealData: replay requires preprocessing_params.slice_normalize.');
        end
        % Require background region for all normalizations (drawn once, reused)
        if ~isfield(preprocessing_params, 'normalize_region') || isempty(preprocessing_params.normalize_region)
            error('preprocessRealData: replay requires preprocessing_params.normalize_region.');
        end
        if flags(strcmp(step_names,'slice_select')) && (~isfield(preprocessing_params,'slices') || isempty(preprocessing_params.slices))
            error('preprocessRealData: replay with slice_select requires preprocessing_params.slices.');
        end
        if flags(strcmp(step_names,'bragg_remove')) && (~isfield(preprocessing_params,'bragg') || ~isfield(preprocessing_params.bragg,'mask2d') || isempty(preprocessing_params.bragg.mask2d))
            error('preprocessRealData: replay with bragg_remove requires preprocessing_params.bragg.mask2d.');
        end
        if flags(strcmp(step_names,'auto_streak')) && (~isfield(preprocessing_params,'autoStreakRemoval_slices') || isempty(preprocessing_params.autoStreakRemoval_slices))
            error('preprocessRealData: replay with auto_streak requires preprocessing_params.autoStreakRemoval_slices.');
        end
        if (flags(strcmp(step_names,'manual_streak')) || flags(strcmp(step_names,'auto_streak')))
            if isfield(preprocessing_params,'streak_params') && ~isempty(preprocessing_params.streak_params)
                % New format: use streak_params
            elseif flags(strcmp(step_names,'auto_streak')) && isfield(preprocessing_params,'autoStreakRemoval_slices') && ~isempty(preprocessing_params.autoStreakRemoval_slices)
                % Backward compat: build streak_params from old fields
                preprocessing_params.streak_params = struct();
                preprocessing_params.streak_params.slices = preprocessing_params.autoStreakRemoval_slices;
                preprocessing_params.streak_params.max_streak_width = 3;
                preprocessing_params.streak_params.mode = 'valley';
                preprocessing_params.streak_params.remove_factor = 1;
                if isfield(preprocessing_params,'autoStreakRemoval_factor1')
                    preprocessing_params.streak_params.remove_factor = preprocessing_params.autoStreakRemoval_factor1;
                end
                preprocessing_params.streak_params.interpolate_factor = 1;
                if isfield(preprocessing_params,'autoStreakRemoval_factor2')
                    preprocessing_params.streak_params.interpolate_factor = preprocessing_params.autoStreakRemoval_factor2;
                end
            else
                error('preprocessRealData: replay with manual_streak or auto_streak requires preprocessing_params.streak_params (or legacy autoStreakRemoval_slices).');
            end
        end
        if flags(strcmp(step_names,'defect_mask'))
            if ~isfield(preprocessing_params,'defect_slice') || ~isfield(preprocessing_params,'num_defect_type')
                error('preprocessRealData: replay with defect_mask requires defect_slice and num_defect_type.');
            end
            if ~isfield(preprocessing_params,'defect_mask') || isempty(preprocessing_params.defect_mask)
                error('preprocessRealData: replay with defect_mask requires preprocessing_params.defect_mask (tg method only).');
            end
        end
        if flags(strcmp(step_names,'heal')) && (~isfield(preprocessing_params,'heal_direction') || isempty(preprocessing_params.heal_direction))
            error('preprocessRealData: replay with heal requires preprocessing_params.heal_direction.');
        end
    else
        % Interactive: defaults for do_* and fresh step flags
        cfgFields = { ...
            "do_bragg_remove", true; "do_crop", true; "do_slice_select", true; ...
            "do_manual_streak", true; "do_auto_streak", true; "do_defect_mask", true; ...
            "do_streak_correct", true; "do_heal", true; "do_directional_plane", false ...
        };
        for i = 1:size(cfgFields,1)
            if ~isfield(params.preprocessing, cfgFields{i,1})
                params.preprocessing.(cfgFields{i,1}) = cfgFields{i,2};
            end
        end
        preprocessing_params = struct();
        preprocessing_params.steps.names = step_names;
        preprocessing_params.steps.flags = false(1, nSteps);
    end

    data_carried = data_original;
    rangetype    = 'dynamic';
    % Use first slice of current stack for normalization (index 1), so valid after slice selection.
    NORM_SLICE = 1;
    preprocessing_params.slice_normalize = 1;
    if ~interactive && isfield(data.real, 'preprocessing') && isfield(data.real.preprocessing, 'normalize_region')
        preprocessing_params.normalize_region = data.real.preprocessing.normalize_region;
    elseif ~isfield(preprocessing_params, 'normalize_region')
        preprocessing_params.normalize_region = [];
    end


    % --- 2.1 Bragg removal ---
    run_step = getRunStep('bragg_remove', interactive, step_names, preprocessing_params, params);
    if run_step
        [data_carried, preprocessing_params] = applyNormalize(data_carried, rangetype, NORM_SLICE, preprocessing_params);
        if interactive
            [data_braggremoved, mask2d, recipe] = removeBragg(data_carried);
            preprocessing_params.bragg = struct('mask2d', mask2d, 'recipe', recipe);
            preprocessing_params.steps.flags(strcmp(step_names,'bragg_remove')) = true;
        else
            opts = struct('mask2d', preprocessing_params.bragg.mask2d, 'recipe', preprocessing_params.bragg.recipe);
            data_braggremoved = removeBragg(data_carried, opts);
        end
        data_carried = data_braggremoved;
        if ~interactive
            preprocessing_params.steps.flags(strcmp(step_names,'bragg_remove')) = true;
        end
    end

    % --- 2.2 Crop ---
    run_step = getRunStep('crop', interactive, step_names, preprocessing_params, params);
    if run_step
        mask = maskSquare(data_carried, 0, 40, 'square');
        data_cropped = gridCropMask(data_carried, mask);
        data_carried = data_cropped;
        preprocessing_params.steps.flags(strcmp(step_names,'crop')) = true;
    else
        data_cropped = data_carried;
    end

    % --- 2.3 Slice selection ---
    run_step = getRunStep('slice_select', interactive, step_names, preprocessing_params, params);
    if run_step
        if interactive
            rangetype = 'dynamic';
            figure;
            d3gridDisplay(data_carried, rangetype);
            preprocessing_params.slices = input('input a list of slices: ');
            params.preprocessing.slices = preprocessing_params.slices;
        else
            preprocessing_params.slices = data.real.preprocessing.slices;
        end
        data_selected = data_carried(:,:,preprocessing_params.slices);
        energy_selected = energy_range(preprocessing_params.slices);
        if interactive, close; end
        data_carried = data_selected;
        preprocessing_params.steps.flags(strcmp(step_names,'slice_select')) = true;
    else
        energy_selected = energy_range;
        preprocessing_params.slices = 1:numel(energy_range);
    end

    % --- 2.4 Streak removal (one run prompt: manual on first slice, then auto on all, repeat/end) ---
    % Replay: preprocessRealData calls wrapper with interactive=false and streak_params → wrapper runs auto only.
    run_streak = getRunStreakStep(interactive, step_names, preprocessing_params, params);
    if run_streak
        if interactive
            opts = struct();
        else
            opts = struct('streak_params', preprocessing_params.streak_params);
        end
        [data_carried, streak_params] = streakRemovalWorkflow(data_carried, interactive, opts);
        preprocessing_params.streak_params = streak_params;
        preprocessing_params.steps.flags(strcmp(step_names,'manual_streak')) = true;
        preprocessing_params.steps.flags(strcmp(step_names,'auto_streak')) = true;
    end

    % --- 2.5 Defect masking ---
    run_step = getRunStep('defect_mask', interactive, step_names, preprocessing_params, params);
    if run_step
        [data_carried, preprocessing_params] = applyNormalize(data_carried, rangetype, NORM_SLICE, preprocessing_params);
        if interactive
            f1 = figure;
            d3gridDisplay(data_carried, 'dynamic');
            preprocessing_params.defect_slice = input('Enter defect slice number: ');
            preprocessing_params.num_defect_type = input('enter how many types of defects to mask: ');
            close(f1);
            params.preprocessing.defect_slice = preprocessing_params.defect_slice;
            params.preprocessing.num_defect_type = preprocessing_params.num_defect_type;
        else
            preprocessing_params.defect_slice = data.real.preprocessing.defect_slice;
            preprocessing_params.num_defect_type = data.real.preprocessing.num_defect_type;
        end

        % Defect masking: tg (Gaussian) only in this wrapper; gw/threshold remain in legacy scripts.
        if interactive
            [data_masked, preprocessing_params.defect_mask, ~, ~] = ...
                gaussianMaskDefects(data_carried, preprocessing_params.defect_slice, preprocessing_params.num_defect_type);
        else
            data_masked = gaussianMaskDefects(data_carried, preprocessing_params.defect_slice, ...
                preprocessing_params.num_defect_type, data.real.preprocessing.defect_mask);
        end
        data_carried = data_masked;
        preprocessing_params.steps.flags(strcmp(step_names,'defect_mask')) = true;
    else
        data_masked = data_carried;
    end

    % --- 2.6a Streak correct ---
    run_step = getRunStep('streak_correct', interactive, step_names, preprocessing_params, params);
    if run_step
        [data_carried, preprocessing_params] = applyNormalize(data_carried, rangetype, NORM_SLICE, preprocessing_params);
        [data_streakremoved, ~] = RemoveStreaks(data_carried, 'Direction', 'vertical');
        data_carried = data_streakremoved;
        preprocessing_params.steps.flags(strcmp(step_names,'streak_correct')) = true;
    end

    % --- 2.6b Heal ---
    run_step = getRunStep('heal', interactive, step_names, preprocessing_params, params);
    if run_step
        [data_carried, preprocessing_params] = applyNormalize(data_carried, rangetype, NORM_SLICE, preprocessing_params);
        if interactive
            preprocessing_params.heal_direction = input('Enter direction to heal (horizontal/vertical/none): ', 's');
            params.preprocessing.heal_direction = preprocessing_params.heal_direction;
        else
            preprocessing_params.heal_direction = data.real.preprocessing.heal_direction;
        end
        data_streakremoved_healed = heal_streaks(data_carried, preprocessing_params.heal_direction);
        data_carried = data_streakremoved_healed;
        preprocessing_params.steps.flags(strcmp(step_names,'heal')) = true;
    end

    % --- 2.6c Directional plane (retired) ---
    run_step = getRunStep('directional_plane', interactive, step_names, preprocessing_params, params);
    if run_step
        [data_carried, preprocessing_params] = applyNormalize(data_carried, rangetype, NORM_SLICE, preprocessing_params);
        preprocessing_params.real_space_direction = 'horizontal';
        [data_plane, ~] = d3plane_directional(data_carried, preprocessing_params.real_space_direction, 'LineWidth', 5);
        data_carried = data_plane;
        preprocessing_params.steps.flags(strcmp(step_names,'directional_plane')) = true;
    end

    % Final normalization
    [Y, preprocessing_params] = applyNormalize(data_carried, rangetype, NORM_SLICE, preprocessing_params);

    data.real.Y = Y;
    data.real.data_cropped = [];
    data.real.data_masked = [];
    if exist('data_cropped','var') && ~isempty(data_cropped)
        data.real.data_cropped = data_cropped;
    end
    if exist('data_masked','var') && ~isempty(data_masked)
        data.real.data_masked = data_masked;
    end
    data.real.energy_selected = energy_selected;
    data.real.preprocessing = preprocessing_params;

    params.preprocessing.steps = preprocessing_params.steps;
    meta.stage = "preprocess";
end

function runStep = getRunStep(stepName, interactive, step_names, preprocessing_params, params)
    idx = strcmp(step_names, stepName);
    if interactive
        doField = ['do_' stepName];
        defaultRun = true;
        if isfield(params.preprocessing, doField)
            defaultRun = params.preprocessing.(doField);
        end
        runStep = askRunOrSkip(stepName, true, defaultRun);
    else
        runStep = preprocessing_params.steps.flags(idx);
    end
end

function runStep = getRunStreakStep(interactive, step_names, preprocessing_params, params)
% Single run prompt for streak removal (manual then auto, repeat/end). Replay = run wrapper with auto only.
    if interactive
        defaultRun = true;
        if isfield(params.preprocessing, 'do_manual_streak')
            defaultRun = params.preprocessing.do_manual_streak;
        end
        if isfield(params.preprocessing, 'do_auto_streak')
            defaultRun = defaultRun || params.preprocessing.do_auto_streak;
        end
        runStep = askRunOrSkip('streak_removal (manual on first slice, then auto on all, repeat/end)', true, defaultRun);
    else
        idx_manual = strcmp(step_names, 'manual_streak');
        idx_auto   = strcmp(step_names, 'auto_streak');
        runStep = preprocessing_params.steps.flags(idx_manual) | preprocessing_params.steps.flags(idx_auto);
    end
end

function runStep = askRunOrSkip(processName, interactive, defaultRun)
    if interactive
        fprintf('Process: %s\n', processName);
        reply = input('Run? (y/skip): ', 's');
        runStep = strcmpi(strtrim(reply), 'y') || strcmpi(strtrim(reply), 'yes');
    else
        runStep = defaultRun;
    end
end

function [data_out, pparams] = applyNormalize(data_in, rangetype, norm_slice, pparams)
% Normalize using shared background region; acquire region (user draws) once if not set.
    if ~isfield(pparams, 'normalize_region') || isempty(pparams.normalize_region)
        [data_out, ~, ~, pos] = normalizeBackgroundToZeroMean3D(data_in, rangetype, norm_slice);
        pparams.normalize_region = pos;
    else
        data_out = normalizeBackgroundToZeroMean3D(data_in, rangetype, norm_slice, pparams.normalize_region);
    end
end
