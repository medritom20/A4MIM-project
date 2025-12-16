clearvars; clc; close all
addpath(genpath('..\'))

dataFile = 's3dkt3m2.mat';              % 's3dkt3m2.mat' | 'bcsstk03.mat'
[A, ~] = load_problem(dataFile);
n = size(A, 1);

nRHS = 16;                                       % 1 | 2 | 4 | 6 | ...
B = randn(n, nRHS);
Xtrue = A \ B;
X0 = zeros(n, nRHS);

% prec.type = 'none';
prec.type = 'ichol';                            % 'ichol' | 'none'
prec.opts.type = 'ict';
prec.opts.droptol = 1e-5;
prec.opts.diagcomp = 1e-2;
applyM = build_preconditioner(A, prec);         % @(u) M^{-1}u

maxIt = 500;
tol = 1e-14;

[X, omegaHist] = PDPBCG(A, B, X0, applyM, Xtrue, maxIt, tol);

plot_convergence(omegaHist, get_legend_label(prec) );
