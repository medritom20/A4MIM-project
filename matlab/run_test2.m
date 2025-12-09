clearvars; clc; close all
dbstop if error
addpath(genpath('.'));

%% EXPERIMENT 2: Reproduce results from Fig. 2 in the paper

matrixName  = 's3dkt3m2.mat';
blockSizes  = [1, 4, 16, 64];
maxIts      = [3500, 1200, 600, 200];
tol         = 1e-10;
rngSeed     = 7;        % Pick a positive integer seed for reproducibility
precCfg = struct('type', 'ichol', 'opts', struct('type', 'ict', 'droptol', 1e-5, 'diagcomp', 1e-2));
    % precCfg = struct('type', 'none');

sessions = repmat(struct(), numel(blockSizes), 1);
for iS = 1:numel(blockSizes)
    sessions(iS).matrixName = matrixName;
    sessions(iS).blockSize = blockSizes(iS);
    sessions(iS).rngSeed = rngSeed + (iS-1);        % Provide session-specific rngSeed for reproducibility
    sessions(iS).prec = precCfg;

    basePolicy = struct('tol', tol, 'maxIt', maxIts(iS));

    % 1st run: DP-BCG
    sessions(iS).runs(1) = struct('method', 'DP-BCG', 'policy', basePolicy);
    % 2nd run: DR-BCG
    sessions(iS).runs(2) = struct('method', 'DR-BCG', 'policy', basePolicy);
    % 3rd run: HS-BCG
    hsPolicy = basePolicy;
    hsPolicy.phiType = 'eye';
    sessions(iS).runs(3) = struct('method', 'HS-BCG', 'policy', hsPolicy);
end

results = driver(sessions);

%% TEMPORARY CODE SEGMENT: LOAD SAVED RESULTS
% clc; close all
% dbstop if error
% load('test2_results.mat')

%% PLOTTING RESULTS
% Prepare plot configuration
styleDefs = method_registry(precCfg);
plotCfg.layout = [2, 2];
plotCfg.title = 'Experiment 2 reproduction: s3dkt3m2 (ICT)';
plotCfg.outputDir = fullfile(thisDir, 'out', 'results');
plotCfg.figureSlug = 'fig2_s3dkt3m2';
plotCfg.saveOutputs = true;

% Generate plots (+SAVE)
plot_session_results(results, styleDefs, plotCfg);
