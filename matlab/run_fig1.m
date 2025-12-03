clearvars; clc; close all
dbstop if error

thisDir = fileparts(mfilename('fullpath'));
addpath(thisDir);
addpath(fullfile(thisDir, 'aux_funs'));
addpath(genpath(fullfile(thisDir, 'methods')));

policy.problemFile = fullfile(thisDir, 'data', 'bcsstk03.mat');
policy.blockSizes  = [1, 2, 4, 6];
policy.methods     = {'HS-BCG', 'DP-BCG', 'DR-BCG'};
policy.maxIt       = 1000;
policy.tol         = 1e-12;
policy.rngSeed     = 42;
policy.prec        = struct('type', 'none');
policy.outputDir   = fullfile(thisDir, 'out', 'results');
policy.figureSlug  = 'fig1_bcsstk03';
policy.figureTitle = 'bcsstk03';
% policy.saveOutputs = true ;

driver(policy);
