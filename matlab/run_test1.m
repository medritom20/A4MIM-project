clearvars; clc; close all

addpath(genpath('.'));

matrixName = 'bcsstk03.mat';
blockSizes = [1, 2, 4, 6];
maxIts = [800, 500, 350, 250];
tol = 1e-10;
rngSeed = 42;

sessions = repmat(struct(), numel(blockSizes), 1);
for iBS = 1:numel(blockSizes)
    sessions(iBS).matrixName = matrixName;
    sessions(iBS).blockSize = blockSizes(iBS);
    sessions(iBS).rngSeed = rngSeed + iBS - 1;
    sessions(iBS).prec = struct('type', 'none');

    basePolicy = struct('tol', tol, 'maxIt', maxIts(iBS));
    % Run 1 - HS-BCG
    hsPolicy = basePolicy;
    hsPolicy.phiType = 'eye';
    sessions(iBS).runs(1) = struct('method', 'HS-BCG', 'policy', hsPolicy);
    % Run 2 - DP-BCG
    sessions(iBS).runs(2) = struct('method', 'DP-BCG', 'policy', basePolicy);
    % Run 3 - DR-BCG
    sessions(iBS).runs(3) = struct('method', 'DR-BCG', 'policy', basePolicy);
end

results = driver(sessions);

styleDefs = method_registry(results.precInfo);
plotCfg.layout = [2, 2];
plotCfg.title = 'Fig. 1 reproduction: bcsstk03';
plotCfg.outputDir = fullfile('.', 'out', 'results');
plotCfg.figureSlug = 'fig1_bcsstk03';
plotCfg.saveOutputs = true;

plot_session_results(results, styleDefs, plotCfg);
