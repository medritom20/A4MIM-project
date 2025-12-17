clearvars; clc; close all
dbstop if error
addpath(genpath('.'));

%% EXPERIMENT #3: DP-BCG vs. BF-BCG (Fig. 3 reproduction)
    info.experimentId = 3;
    info.rngSeed      = 50;

    solver.matrixName = 's3dkt3m2.mat';
    solver.blockSizes = [16, 64];
    solver.maxIts     = [500, 200];
    solver.tol        = 1e-15;
    solver.precCfg    = struct('type', 'ichol', 'opts', struct('type', 'ict', ...
                         'droptol', 1e-5, 'diagcomp', 1e-2));
    solver.svdTols     = [  1e-3, 1e-4, 1e-5, 1e-6, ...
                            1e-7, 1e-8, 1e-9, 1e-10];
    solver.bfStyles   = {   ':', '--', '-.', '-', ...
                            ':', '--', '-.', '-'};

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
clc
    % N shades: dark red -> orange -> yellow (no black/white)
    N = 8;
    
    M  = 256;                 % resolution of base map
    cm = hot(M);              % [black -> red -> yellow -> white]
    
    % skip black-ish part &  avoid near-white end
    idx = round(linspace(round(0.15*M), round(0.65*M), numel(solver.svdTols))); 
    
    colors = cm(idx, :);      % Nx3 RGB in [0,1]
    
    solver.colors = colors ;

styleDefs = method_registry(solver);
plotCfg.layout = [2, 1];
plotCfg.title = 'DP-BCG vs. BF-BCG (s3dkt3m2, ICT)';
plotCfg.outputDir = fullfile('.\out', 'figs','EXP03');
plotCfg.figureSlug = sprintf('fig3_seed%03i',info.rngSeed);
plotCfg.saveOutputs = true;
plotCfg.experimentId = info.experimentId;


plot_session_results(results, styleDefs, plotCfg);


