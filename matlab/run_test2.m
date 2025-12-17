clearvars; clc; close all
dbstop if error
addpath(genpath('.'));

%% EXPERIMENT #2: Reproduce results from Fig. 2 in the paper
    % Mirror run_test1 structure for consistency
    info.experimentId = 2;
    info.rngSeed      = 10;

    solver.matrixName = 's3dkt3m2.mat';
    solver.blockSizes = [1, 4, 16, 64];
    solver.maxIts     = [3500, 1200, 600, 200];
    solver.tol        = 1e-10;
    solver.precCfg    = struct('type', 'ichol', 'opts', struct('type', 'ict', ...
                        'droptol', 1e-5, 'diagcomp', 1e-2));

    % SETUP SESSIONS
    sessions = repmat(struct(), numel(solver.blockSizes), 1);
    for iS = 1:numel(solver.blockSizes)
        sessions(iS).matrixName = solver.matrixName;
        sessions(iS).blockSize  = solver.blockSizes(iS);
        sessions(iS).rngSeed    = info.rngSeed + (iS-1);
        sessions(iS).prec       = solver.precCfg;

        basePolicy = struct('tol', solver.tol, 'maxIt', solver.maxIts(iS));

        % 1st run: HS-BCG
        hsPolicy = basePolicy;
        hsPolicy.phiType = 'eye';
        sessions(iS).runs(1) = struct('method', 'HS-BCG', 'policy', hsPolicy);
        % 2nd run: DP-BCG
        sessions(iS).runs(2) = struct('method', 'DP-BCG', 'policy', basePolicy);
        % 3rd run: DR-BCG
        sessions(iS).runs(3) = struct('method', 'DR-BCG', 'policy', basePolicy);
    end

results = driver(sessions);

%% SAVE RESULTS TO TABLE
    outputTableDir = ensure_dir( fullfile('.\out', 'results') );
    save_results( solver.matrixName, info.rngSeed, solver.precCfg, results, sessions, info.experimentId);

%% PLOTTING RESULTS
% Prepare plot configuration
styleDefs = method_registry(solver);
plotCfg.layout = [2, 2];
plotCfg.title = 'Experiment 2 reproduction: s3dkt3m2 (ICT)';
plotCfg.outputDir = fullfile('.\out', 'results');
plotCfg.figureSlug = 'fig2_s3dkt3m2';
plotCfg.saveOutputs = true;
plotCfg.experimentId = info.experimentId;

% Generate plots (+SAVE)
plot_session_results(results, styleDefs, plotCfg);
