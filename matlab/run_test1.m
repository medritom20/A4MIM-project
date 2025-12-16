clearvars; clc; close all
dbstop if error
addpath(genpath('.'));

%% EXPERIMENT #1: Reproduce results from Fig. 1 in the paper
    % Use structures to keep the code organized
    info.experimentId   = 1 ; % 1st experiment, do not change!
    info.rngSeed        = 1;  % Provide session-specific rngSeed for reproducibility

    solver.matrixName  = 'bcsstk03.mat'; 
    solver.blockSizes  = [1, 2, 4, 6];
    solver.maxIts      = [1000, 600, 400, 200];
    solver.tol         = 1e-15; 
    solver.precCfg     = struct('type', 'none');

    % SETUP SESSIONS
    sessions = repmat(struct(), numel(solver.blockSizes), 1);
    for iS = 1:numel(solver.blockSizes)
        sessions(iS).matrixName = solver.matrixName;
        sessions(iS).blockSize  = solver.blockSizes(iS);
        sessions(iS).rngSeed    = solver.rngSeed + (iS-1);    % set unique rngSeed for each run to avoid identical solutions  
        sessions(iS).prec       = solver.precCfg; 

        basePolicy = struct('tol', solver.tol, 'maxIt', solver.maxIts(iS));
        % Run 1 - HS-BCG
        hsPolicy = basePolicy;
        hsPolicy.phiType = 'eye';
        sessions(iS).runs(1) = struct('method', 'HS-BCG', 'policy', hsPolicy);
        % Run 2 - DP-BCG
        sessions(iS).runs(2) = struct('method', 'DP-BCG', 'policy', basePolicy);
        % Run 3 - DR-BCG
        sessions(iS).runs(3) = struct('method', 'DR-BCG', 'policy', basePolicy);
    end

results = driver(sessions);

%% SAVE RESULTS TO TABLE
    outputTableDir = ensure_dir( fullfile('.\out', 'results') );
    save_results( solver.matrixName, info.rngSeed, solver.precCfg, results, sessions, experimentId);

%% PLOTTING RESULTS
    % Prepare plot configuration
    styleDefs = method_registry(solver);
    plotCfg.layout = [2, 2];
    plotCfg.title = 'Experiment 1 reproduction: bcsstk03';
    plotCfg.outputDir = fullfile('.\out', 'results');
    plotCfg.figureSlug = 'fig1_bcsstk03';
    plotCfg.saveOutputs = true;
    plotCfg.experimentId = info.experimentId;

    % Generate plots (+SAVE)
    plot_session_results(results, styleDefs, plotCfg);
