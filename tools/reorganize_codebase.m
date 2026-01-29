function reorganize_codebase(varargin)
%REORGANIZE_CODEBASE  Apply the repo reorg described in Reorg_plan.md
%
% Usage (from repo root):
%   reorganize_codebase('dryRun', true)   % prints actions only (default)
%   reorganize_codebase('dryRun', false)  % actually moves/edits files
%
% This function is intentionally conservative:
% - It refuses to overwrite existing destinations.
% - It prints every action it would take.
%
% It implements Steps 1–9 from `Reorg_plan.md`:
% - Create new directory skeleton
% - Move third-party libraries into lib/
% - Move helpers/solvers/algorithms/wrappers/etc. into src/
% - Move tests/archive/tools
% - Apply the small set of path/text fixes (init_sbd, startup, solvers, phasespace)

    p = inputParser;
    p.addParameter('dryRun', true, @(x) islogical(x) && isscalar(x));
    p.parse(varargin{:});
    dryRun = p.Results.dryRun;

    thisFile = mfilename('fullpath');              % .../tools/reorganize_codebase.m
    repoRoot = fileparts(fileparts(thisFile));     % .../ (repo root)

    if ~exist(fullfile(repoRoot, 'Reorg_plan.md'), 'file')
        error('Could not locate repo root (Reorg_plan.md missing). Run from within this repository.');
    end

    fprintf('\n== MT-SBD-STM reorganization ==\n');
    fprintf('Repo root: %s\n', repoRoot);
    fprintf('Mode: %s\n\n', ternary(dryRun, 'DRY RUN (no changes)', 'APPLY (will change files)'));

    %% Step 1: Create directory skeleton
    dirs = {
        'src/algorithms'
        'src/workflow'
        'src/solvers'
        'src/helpers'
        'src/initialization'
        'src/metrics'
        'src/generation'
        'src/preprocessing'
        'src/io'
        'src/visualization'
        'src/utils'
        'tests'
        'lib'
        'lib/imshow3D'
        'lib/colormap'
        'lib/mat2im'
        'archive'
        'tools'
    };
    for i = 1:numel(dirs)
        ensureDir(fullfile(repoRoot, dirs{i}), dryRun);
    end

    %% Step 2: Move third-party libraries to lib/
    % We move *contents* into lib/* so it works even if lib/* exists.
    moveFolderContentsIfExists(fullfile(repoRoot, 'Dong_func', 'imshow3D.m'), fullfile(repoRoot, 'lib', 'imshow3D'), dryRun);
    removeDirIfEmpty(fullfile(repoRoot, 'Dong_func', 'imshow3D.m'), dryRun);

    moveFolderContentsIfExists(fullfile(repoRoot, 'colormap'), fullfile(repoRoot, 'lib', 'colormap'), dryRun);
    removeDirIfEmpty(fullfile(repoRoot, 'colormap'), dryRun);

    moveFolderContentsIfExists(fullfile(repoRoot, 'utils', 'mat2im'), fullfile(repoRoot, 'lib', 'mat2im'), dryRun);
    removeDirIfEmpty(fullfile(repoRoot, 'utils', 'mat2im'), dryRun);

    %% Step 3: Move core helpers to src/helpers/
    helperFiles = {
        'cconvfft2.m'
        'convfft2.m'
        'convfft3.m'
        'extend.m'
        'H_function.m'
        'huber.m'
        'Hxx_function.m'
        'proj2oblique.m'
        'seq_crosscorr_regularizer.m'
        'shrink.m'
    };
    for i = 1:numel(helperFiles)
        moveFileIfExists( ...
            fullfile(repoRoot, 'core', 'helpers', helperFiles{i}), ...
            fullfile(repoRoot, 'src', 'helpers', helperFiles{i}), ...
            dryRun);
    end

    %% Step 4: Move solvers to src/solvers/ and tests/
    solverFiles = {
        'Asolve_Manopt.m'
        'Asolve_Manopt_ALL.m'
        'Asolve_Manopt_multikernel.m'
        'Asolve_Manopt_parallel.m'
        'Asolve_Manopt_tunable.m'
        'trustregions.m'
        'Xsolve_FISTA.m'
        'Xsolve_FISTA_ALL.m'
        'Xsolve_FISTA_multikernel.m'
        'Xsolve_FISTA_parallel.m'
        'Xsolve_FISTA_tunable.m'
        'Xsolve_pdNCG.m'
    };
    for i = 1:numel(solverFiles)
        moveFileIfExists( ...
            fullfile(repoRoot, 'core', solverFiles{i}), ...
            fullfile(repoRoot, 'src', 'solvers', solverFiles{i}), ...
            dryRun);
    end

    % Apply the required "depth changed" path fixes in the moved solver files.
    fixSolverConfigPaths(fullfile(repoRoot, 'src', 'solvers'), dryRun);

    % Move solver tests out of src
    moveFileIfExists(fullfile(repoRoot, 'core', 'Asolve_Manopt_tunable_test.m'), fullfile(repoRoot, 'tests', 'Asolve_Manopt_tunable_test.m'), dryRun);
    moveFileIfExists(fullfile(repoRoot, 'core', 'Xsolve_FISTA_test.m'),         fullfile(repoRoot, 'tests', 'Xsolve_FISTA_test.m'),         dryRun);

    % Centralize tunable configs if they exist.
    % If a different config already exists at the destination, don't overwrite:
    % move the incoming one into config/legacy/ with a suffix instead.
    moveFileIfExists(fullfile(repoRoot, 'examples', 'Asolve_config_tunable.mat'), fullfile(repoRoot, 'config', 'Asolve_config_tunable.mat'), dryRun, ...
        'onConflict', 'rename', 'conflictTag', 'from_examples');
    moveFileIfExists(fullfile(repoRoot, 'examples', 'Xsolve_config_tunable.mat'), fullfile(repoRoot, 'config', 'Xsolve_config_tunable.mat'), dryRun, ...
        'onConflict', 'rename', 'conflictTag', 'from_examples');
    moveFileIfExists(fullfile(repoRoot, 'Asolve_config_tunable.mat'),            fullfile(repoRoot, 'config', 'Asolve_config_tunable.mat'), dryRun, ...
        'onConflict', 'rename', 'conflictTag', 'from_root');
    moveFileIfExists(fullfile(repoRoot, 'Xsolve_config_tunable.mat'),            fullfile(repoRoot, 'config', 'Xsolve_config_tunable.mat'), dryRun, ...
        'onConflict', 'rename', 'conflictTag', 'from_root');

    %% Step 5: Move algorithm entrypoints to src/algorithms/
    algoFiles = {
        'MT_SBD.m'
        'SBD.m'
        'MTSBD_synthetic.m'
        'MTSBD_synthetic_Xregulated.m'
        'MTSBD_synthetic_all_slice.m'
        'MTSBD_synthetic_Xregulated_all_slices.m'
        'MTSBD_all_slice.m'
    };
    for i = 1:numel(algoFiles)
        moveFileIfExists( ...
            fullfile(repoRoot, algoFiles{i}), ...
            fullfile(repoRoot, 'src', 'algorithms', algoFiles{i}), ...
            dryRun);
    end

    %% Step 6: Decompose Dong_func/ into src/ modules
    % Subfolders with explicit targets
    moveFolderContentsIfExists(fullfile(repoRoot, 'Dong_func', 'wrapper'), fullfile(repoRoot, 'src', 'workflow'), dryRun);
    moveFolderContentsIfExists(fullfile(repoRoot, 'Dong_func', 'basic'),   fullfile(repoRoot, 'src', 'utils'),    dryRun);
    moveFolderContentsIfExists(fullfile(repoRoot, 'Dong_func', 'data_preprocessing'), fullfile(repoRoot, 'src', 'preprocessing'), dryRun);

    % Individual files (heuristic routing; conservative, but covers all remaining Dong_func/*.m)
    routeDongFuncRootFiles(repoRoot, dryRun);

    %% Step 7: Move remaining root/legacy files and clean up folders
    % Root tests -> tests/
    rootTestFiles = {
        'SBD_test.m'
        'SBD_test_multi.m'
        'SBD_test_multi_overlap.m'
        'SBD_test_multi_parallel.m'
        'SBD_test_multi_demixing_attempt.m'
        'SBD_with_proximity_script.m'
    };
    for i = 1:numel(rootTestFiles)
        moveFileIfExists(fullfile(repoRoot, rootTestFiles{i}), fullfile(repoRoot, 'tests', rootTestFiles{i}), dryRun);
    end

    moveFileIfExists(fullfile(repoRoot, 'plot_activations.m'),   fullfile(repoRoot, 'src', 'visualization', 'plot_activations.m'), dryRun);
    moveFileIfExists(fullfile(repoRoot, 'streak_correction.m'),  fullfile(repoRoot, 'src', 'preprocessing', 'streak_correction.m'), dryRun);
    moveFileIfExists(fullfile(repoRoot, 'function_path_analyzer.m'), fullfile(repoRoot, 'tools', 'function_path_analyzer.m'), dryRun);

    moveFileIfExists(fullfile(repoRoot, 'plotting', 'streaknoise.m'), fullfile(repoRoot, 'src', 'visualization', 'streaknoise.m'), dryRun);
    % historical/* -> archive/ (archive already exists from skeleton)
    moveFolderContentsIfExists(fullfile(repoRoot, 'historical'), fullfile(repoRoot, 'archive'), dryRun);
    removeDirIfEmpty(fullfile(repoRoot, 'historical'), dryRun);

    % Legacy utils/*.m -> src/utils (mat2im already moved)
    % Only move MATLAB source files; leave any stray .mat in place.
    moveFolderContentsIfExists(fullfile(repoRoot, 'utils'), fullfile(repoRoot, 'src', 'utils'), dryRun, ...
        'ignoreDirs', {'mat2im'}, 'fileExtWhitelist', {'.m'});

    %% Step 8: Update init_sbd.m and startup.m
    updateInitSbd(fullfile(repoRoot, 'init_sbd.m'), dryRun);
    updateStartup(fullfile(repoRoot, 'startup.m'), dryRun);

    %% Step 9: Fix init_sbd references in examples
    fixRunInitSbd(fullfile(repoRoot, 'examples', 'phase_space', 'phasespace_dataset.m'), dryRun);

    %% Step 7b: Remove now-empty old directories (best-effort)
    removeDirIfEmpty(fullfile(repoRoot, 'core', 'helpers'), dryRun);
    removeDirIfEmpty(fullfile(repoRoot, 'core'), dryRun);
    removeDirIfEmpty(fullfile(repoRoot, 'colormap'), dryRun);
    removeDirIfEmpty(fullfile(repoRoot, 'plotting'), dryRun);
    removeDirIfEmpty(fullfile(repoRoot, 'Dong_func'), dryRun);
    removeDirIfEmpty(fullfile(repoRoot, 'utils'), dryRun);

    fprintf('\n== Done ==\n');
    if dryRun
        fprintf('No changes were made (dryRun=true).\n');
        fprintf('Next: run reorganize_codebase(''dryRun'', false) to apply.\n');
    else
        fprintf('Changes applied.\n');
        fprintf('Next: in MATLAB, run init_sbd; then `which MT_SBD` / `which convfft2` to verify paths.\n');
    end
end

%% --- helpers ------------------------------------------------------------

function out = ternary(cond, a, b)
    if cond; out = a; else; out = b; end
end

function ensureDir(d, dryRun)
    if exist(d, 'dir'); return; end
    if dryRun
        fprintf('MKDIR %s\n', d);
    else
        mkdir(d);
    end
end

function moveFileIfExists(src, dst, dryRun, varargin)
    ip = inputParser;
    ip.addParameter('onConflict', 'error', @(s) ischar(s) || isstring(s)); % 'error' | 'rename'
    ip.addParameter('conflictTag', '', @(s) ischar(s) || isstring(s));
    ip.parse(varargin{:});
    onConflict = char(ip.Results.onConflict);
    conflictTag = char(ip.Results.conflictTag);

    if exist(src, 'file') ~= 2
        return;
    end
    if exist(dst, 'file') == 2
        % If the destination already exists, only allow it if contents match
        % (e.g., rerunning after a partial move).
        if filesAreIdenticalText(src, dst)
            if dryRun
                fprintf('SKIP  %s (already present at destination)\n', src);
            else
                delete(src);
            end
            return;
        end
        % Different contents.
        if dryRun
            fprintf('CONFLICT %s\n  ->  %s\n', src, dst);
            fprintf('  (destination exists with different contents; would not overwrite)\n');
            return;
        end

        if strcmpi(onConflict, 'rename')
            alt = conflictDestination(dst, conflictTag);
            ensureDir(fileparts(alt), dryRun);
            fprintf('CONFLICT-RENAME %s\n  ->  %s\n', src, alt);
            [ok, msg] = movefile(src, alt);
            if ~ok; error('%s', msg); end
            return;
        end

        error('Refusing to overwrite existing file (different contents): %s', dst);
    end
    ensureDir(fileparts(dst), dryRun);
    if dryRun
        fprintf('MOVE  %s\n  ->  %s\n', src, dst);
    else
        [ok, msg] = movefile(src, dst);
        if ~ok; error('%s', msg); end
    end
end

function moveDirIfExists(src, dst, dryRun)
    if exist(src, 'dir') ~= 7
        return;
    end
    if exist(dst, 'dir') == 7
        error('Refusing to overwrite existing directory: %s', dst);
    end
    ensureDir(fileparts(dst), dryRun);
    if dryRun
        fprintf('MOVEDIR %s\n    -> %s\n', src, dst);
    else
        [ok, msg] = movefile(src, dst);
        if ~ok; error('%s', msg); end
    end
end

function removeDirIfEmpty(d, dryRun)
    if exist(d, 'dir') ~= 7
        return;
    end
    listing = dir(d);
    listing = listing(~ismember({listing.name}, {'.','..'}));
    if ~isempty(listing)
        return; % not empty
    end
    if dryRun
        fprintf('RMDIR %s\n', d);
    else
        rmdir(d);
    end
end

function tf = filesAreIdenticalText(a, b)
    % Best-effort equality for .m files (text). If either read fails, assume different.
    try
        ta = fileread(a);
        tb = fileread(b);
        tf = strcmp(ta, tb);
    catch
        tf = false;
    end
end

function alt = conflictDestination(dst, tag)
    % Place conflicting files under config/legacy/ next to their intended destination.
    [d, base, ext] = fileparts(dst);
    legacyDir = fullfile(d, 'legacy');
    if isempty(tag)
        alt = fullfile(legacyDir, [base '__conflict' ext]);
        return;
    end
    alt = fullfile(legacyDir, [base '__' tag ext]);
end

function moveFolderContentsIfExists(srcDir, dstDir, dryRun, varargin)
    ip = inputParser;
    ip.addParameter('ignoreDirs', {}, @(c) iscell(c));
    ip.addParameter('fileExtWhitelist', {}, @(c) iscell(c));
    ip.parse(varargin{:});
    ignoreDirs = ip.Results.ignoreDirs;
    fileExtWhitelist = ip.Results.fileExtWhitelist;

    if exist(srcDir, 'dir') ~= 7
        return;
    end
    ensureDir(dstDir, dryRun);

    listing = dir(srcDir);
    listing = listing(~ismember({listing.name}, {'.','..'}));
    for i = 1:numel(listing)
        name = listing(i).name;
        src = fullfile(srcDir, name);
        dst = fullfile(dstDir, name);

        if listing(i).isdir
            if any(strcmp(name, ignoreDirs))
                continue;
            end
            moveDirIfExists(src, dst, dryRun);
        else
            if ~isempty(fileExtWhitelist)
                [~, ~, ext] = fileparts(name);
                if ~any(strcmpi(ext, fileExtWhitelist))
                    continue;
                end
            end
            moveFileIfExists(src, dst, dryRun);
        end
    end
end

function routeDongFuncRootFiles(repoRoot, dryRun)
    dongRoot = fullfile(repoRoot, 'Dong_func');
    if exist(dongRoot, 'dir') ~= 7
        return;
    end

    listing = dir(fullfile(dongRoot, '*.m'));
    for i = 1:numel(listing)
        name = listing(i).name;
        src = fullfile(dongRoot, name);

        % Skip items already handled explicitly
        if strcmp(name, 'imshow3D.m'); continue; end

        dst = dongFuncDestination(repoRoot, name);
        moveFileIfExists(src, dst, dryRun);
    end
end

function dst = dongFuncDestination(repoRoot, filename)
    % Heuristic routing by filename (best-effort; organization-only).
    n = filename;

    if startsWith(n, 'initialize_') || strcmp(n, 'initialize_parameter_space.m') || strcmp(n, 'selectKernelsInteractive.m')
        dst = fullfile(repoRoot, 'src', 'initialization', n);
        return;
    end

    if startsWith(n, 'generate') || startsWith(n, 'properGen') || startsWith(n, 'sample_parameter_space') || strcmp(n, 'padWithSelectedAreaNoise.m')
        dst = fullfile(repoRoot, 'src', 'generation', n);
        return;
    end

    if contains(n, 'load') || contains(n, 'save') || contains(n, 'parallel') || strcmp(n, 'adaptParallelResults.m') || strcmp(n, 'load3dsall.m')
        dst = fullfile(repoRoot, 'src', 'io', n);
        return;
    end

    if startsWith(n, 'plot_') || startsWith(n, 'visualize') || startsWith(n, 'd3') || contains(n, 'Display') || contains(n, 'imshow')
        dst = fullfile(repoRoot, 'src', 'visualization', n);
        return;
    end

    if contains(n, 'Streak') || startsWith(n, 'remove') || startsWith(n, 'mask') || contains(n, 'normalize') || contains(n, 'threshold') || contains(n, 'windowToKernel') || contains(n, 'gauss') || contains(n, 'filter')
        dst = fullfile(repoRoot, 'src', 'preprocessing', n);
        return;
    end

    if startsWith(n, 'compute') || startsWith(n, 'evaluate') || startsWith(n, 'compare') || contains(n, 'metric') || contains(n, 'Quality') || contains(n, 'SNR') || contains(n, 'Similarity') || contains(n, 'overlap')
        dst = fullfile(repoRoot, 'src', 'metrics', n);
        return;
    end

    dst = fullfile(repoRoot, 'src', 'utils', n);
end

function updateInitSbd(initPath, dryRun)
    if exist(initPath, 'file') ~= 2
        return;
    end

    txt = fileread(initPath);

    % Replace the subdirectory addpath loop with src/config/lib
    pat = sprintf("for d = {'core', 'utils', 'config'}\n        addpath(genpath([fp d{1}]));\n    end");
    repl = sprintf("addpath(genpath([fp 'src']));\n    addpath(genpath([fp 'config']));\n    addpath(genpath([fp 'lib']));");

    if contains(txt, pat)
        newTxt = strrep(txt, pat, repl);
    else
        % Fallback: do a softer rewrite of the whole "Add subdirectories" block.
        nl = sprintf('\n');
        newBlockLines = {
            '    % Add subdirectories to path'
            '    fp = [fileparts(mfilename(''fullpath'')) ''/''];'
            '    addpath(fp);'
            '    addpath(genpath([fp ''src'']));      % All internal source code'
            '    addpath(genpath([fp ''config'']));   % Configuration'
            '    addpath(genpath([fp ''lib'']));      % Third-party libraries'
            ''
        };
        newBlock = strjoin(newBlockLines, nl);
        newTxt = regexprep(txt, "(\s*% Add subdirectories to path\s*[\s\S]*?)(\s*% Apply default config settings\s*)", [newBlock '$2'], 'once');
    end

    writeFileIfChanged(initPath, txt, newTxt, dryRun, 'UPDATE init_sbd.m');
end

function fixSolverConfigPaths(srcSolversDir, dryRun)
    % Implements the "Critical path fixes" from Reorg_plan.md Step 4.
    if exist(srcSolversDir, 'dir') ~= 7
        return;
    end

    files = {
        'Asolve_Manopt.m'
        'Asolve_Manopt_ALL.m'
        'Asolve_Manopt_multikernel.m'
        'Asolve_Manopt_tunable.m'
        'Xsolve_FISTA.m'
        'Xsolve_FISTA_tunable.m'
        'Xsolve_FISTA_ALL.m'
        'Xsolve_FISTA_multikernel.m'
        'Xsolve_FISTA_parallel.m'
        'Xsolve_pdNCG.m'
    };

    for i = 1:numel(files)
        pth = fullfile(srcSolversDir, files{i});
        if exist(pth, 'file') ~= 2
            continue;
        end

        old = fileread(pth);
        new = old;

        % Config path depth changes (core -> src/solvers is one level deeper)
        new = strrep(new, "/../config/", "/../../config/");
        new = strrep(new, "/../examples/Asolve_config_tunable.mat", "/../../config/Asolve_config_tunable.mat");
        new = strrep(new, "/../examples/Xsolve_config_tunable.mat", "/../../config/Xsolve_config_tunable.mat");

        % Remove now-unnecessary "helpers" addpath lines in Xsolve files
        new = regexprep(new, "^\s*addpath\(\[fpath\s*'\s*/helpers'\s*\]\);\s*\r?\n", "", 'lineanchors');

        writeFileIfChanged(pth, old, new, dryRun, 'FIX solver config paths');
    end
end

function updateStartup(startupPath, dryRun)
    if exist(startupPath, 'file') ~= 2
        return;
    end

    txt = fileread(startupPath);
    newTxt = txt;

    % Replace old colormap path with lib/colormap (or lib in general).
    newTxt = regexprep(newTxt, ...
        "addpath\(genpath\(fullfile\(pwd,\s*'colormap'\)\)\);\s*", ...
        "addpath(genpath(fullfile(pwd, 'lib')));  % libs (includes lib/colormap)\n", ...
        'once');

    writeFileIfChanged(startupPath, txt, newTxt, dryRun, 'UPDATE startup.m');
end

function fixRunInitSbd(pathToFile, dryRun)
    if exist(pathToFile, 'file') ~= 2
        return;
    end
    txt = fileread(pathToFile);
    newTxt = strrep(txt, "run('../init_sbd');", "run('../../init_sbd');");
    writeFileIfChanged(pathToFile, txt, newTxt, dryRun, 'FIX phasespace_dataset init_sbd path');
end

function writeFileIfChanged(pathToFile, oldTxt, newTxt, dryRun, label)
    if strcmp(oldTxt, newTxt)
        return;
    end
    if dryRun
        fprintf('%s %s (content change)\n', label, pathToFile);
        return;
    end
    fid = fopen(pathToFile, 'w');
    if fid < 0
        error('Failed to open for writing: %s', pathToFile);
    end
    cleaner = onCleanup(@() fclose(fid));
    fwrite(fid, newTxt);
end

