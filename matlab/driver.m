clearvars; clc; close all
addpath(genpath('.\'))

dataFile = './data/bcsstk03.mat';              % './data/s3dkt3m2.mat' | './data/bcsstk03.mat'
[A, ~] = load_problem(dataFile);
n = size(A, 1);

nRHS = 4;                                       % 1 | 2 | 4 | 6 | ...
Xtrue = randn(n, nRHS);
B = A * Xtrue;
X0 = zeros(n, nRHS);

prec.type = 'ichol';                            % 'ichol' | 'none'
prec.opts.type = 'ict';
prec.opts.droptol = 1e-5;
prec.opts.diagcomp = 1e-2;
applyM = build_preconditioner(A, prec);         % @(u) M^{-1}u

maxIt = 100;
tol = 1e-8;

[X, omegaHist] = PDPBCG(A, B, X0, applyM, Xtrue, maxIt, tol);

plot_convergence(omegaHist, get_legend_label(prec) );

