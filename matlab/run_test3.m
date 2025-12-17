clearvars; clc; close all
dbstop if error
addpath(genpath('.'));

%% EXPERIMENT #3: DP-BCG vs. BF-BCG (Fig. 3 reproduction)
    info.experimentId = 3;
    info.rngSeed      = 10;

    solver.matrixName = 's3dkt3m2.mat';
    solver.blockSizes = [16, 64];
    solver.maxIts     = [500, 200];
    solver.tol        = 1e-15;
    solver.precCfg    = struct('type', 'ichol', 'opts', struct('type', 'ict', ...
                         'droptol', 1e-5, 'diagcomp', 1e-2));
    solver.svdTols     = [1e-7, 1e-8, 1e-9, 1e-10];
    solver.bfStyles   = {':', '--', '-.', '-'};

    sessions = repmat(struct(), numel(solver.blockSizes), 1);
    for iS = 1:numel(solver.blockSizes)
        iRun = 0 ;
        sessions(iS).matrixName = solver.matrixName;
        sessions(iS).blockSize  = solver.blockSizes(iS);
        sessions(iS).rngSeed    = info.rngSeed + (iS-1);
        sessions(iS).prec       = solver.precCfg;

        basePolicy = struct('tol', solver.tol, 'maxIt', solver.maxIts(iS));

        for iSVDTol = 1:numel(solver.svdTols)
            iRun = iRun+1 ;
            % 2nd-5th runs: PBF-BCG
            bfPolicy = basePolicy ;
            bfPolicy.svdTol = solver.svdTols(iSVDTol) ;
            
            sessions(iS).runs(iRun) = struct('method', 'BF-BCG', 'policy', bfPolicy);
        end
        % 1st run: PDP-BCG
        dpPolicy = basePolicy ;
        iRun = iRun+1;
        sessions(iS).runs(iRun) = struct('method', 'DP-BCG', 'policy', dpPolicy);
    end

results = driver(sessions);

%% SAVE RESULTS TO TABLE
    outputTableDir = ensure_dir( fullfile('.\out', 'results') );
    save_results( solver.matrixName, info.rngSeed, solver.precCfg, results, sessions, info.experimentId);

%% PLOTTING RESULTS
styleDefs = method_registry(solver);
plotCfg.layout = [2, 1];
plotCfg.title = 'DP-BCG vs. BF-BCG (s3dkt3m2, ICT)';
plotCfg.outputDir = fullfile('.\out', 'results');
plotCfg.figureSlug = sprintf('fig3_seed%03i',info.rngSeed);
plotCfg.saveOutputs = true;
plotCfg.experimentId = info.experimentId;


plot_session_results(results, styleDefs, plotCfg);
