clearvars; clc; close all
dbstop if error

addpath(genpath('.\'));

%% SOLVER SETUP + POLICY
policyDefault.outputDir   = fullfile(thisDir, 'out', 'results');
policyDefault.figureSlug  = 'fig2_s3dkt3m2';
policyDefault.figureTitle = '';

solverSetup.problemFile = '.\data\s3dkt3m2.mat';
solverSetup.blockSizes  = [1,4, 16, 64];
solverSetup.maxIts      = [3500, 1200, 600, 200];
solverSetup.rngSeed     = 7;

policyDefault.tol         = 1e-12;
policyDefault.prec.type   = 'ichol';
policyDefault.prec.opts   = struct( ...
    'type',     'ict', ...
    'droptol',  1e-5, ...
    'diagcomp', 1e-2);

%% SET UP ORCHESTRATOR
for m = blockSizes
    policy = policyDefault ;
    policy.maxIt = solverSetup.maxIts(m) ;
        sess(m).orchestrator(1) = struct( 'method','HS-BCG', 'policy', policy) ;
        sess(m).orchestrator(2) = struct( 'method','DP-BCG', 'policy', policy) ;
        sess(m).orchestrator(3) = struct( 'method','DR-BCG', 'policy', policy) ;
end

%% CALL METHODS' ORCHESTRATOR AND SOLVE SYSTEMS
driver(sess);


%% POSTPROCESSING, PLOTS, SAVE RESULTS

% tady uděláme nějaký cyklus přes sessions, vytvoříme 2x2 subplot a v cyklu
% zavoláme plotovací funkci (uloženou do samostatného souboru), která
% vykreslí do fig pomocí hold on výsledky každé jednotlivé metody a přidá
% do legendy příslušný popisek; toto nám umožní tuto funkci recyklovat i
% pro další ploty (např. v run_fig3.m atd.) a kód bude mnohem čistší. 

someAuxiliaryStruct.outputDir   = fullfile(thisDir, 'out', 'results');
someAuxiliaryStruct.figureSlug  = 'fig1_bcsstk03';
someAuxiliaryStruct.figureTitle = 'bcsstk03';
% someAuxiliaryStruct.saveOutputs = true ;









