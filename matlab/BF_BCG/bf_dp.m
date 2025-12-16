%%%%%%%%%%%%%%%%
% directory management
%%%%%%%%%%%%%%%%

srcDir = fullfile(pwd, "..");
bfDir = fullfile(srcDir, "BF_BCG");
pdpDir = fullfile(srcDir, "PDP_BCG");
dataDir = fullfile(srcDir, "data", "s3dkt3m2.mat");
outDir = fullfile(srcDir, "out", "results", "dp_bf");
addpath(srcDir);
addpath(bfDir);
addpath(pdpDir);

% loading the data
A = load(dataDir).Problem.A;
n = length(A);

%%%%%%%%%%%%%%%%
% preconditioner options
%%%%%%%%%%%%%%%%

opts.type = 'ict';        % threshold dropping
opts.droptol = 1e-5;      % drop tolerance
opts.diagcomp = 1e-2;     % global diagonal shift

L = ichol(A, opts);       
applyM = @(R) L' \ (L \ R);


% iterations, tolerances
maxit = 500;
tol = 1e-8;
thrs = [1e-7, 1e-8, 1e-9, 1e-10];
linespecs = ["r:", "r-.", "r--", "r-"];

%%%%%%%%%%%%%%%%
%m = 16
%%%%%%%%%%%%%%%%

m1 = 16;
X_ex1 = rand(n,m1);
b1 = A*X_ex1;
X_zero1 = zeros(n,m1);

fig16 = figure;
ax16 = axes(fig16);
hold(ax16, 'on');

for i=1:length(thrs)
	thr = thrs(i);
	linespec = linespecs(i);
	[~, ~, ~, omegahist] = bf_bcg(A, b1, applyM, X_ex1, maxit, tol, thr);
	semilogy(ax16, omegahist, linespec);
end

[~, omegahist] = PDPBCG(A, b1, X_zero1, applyM, X_ex1, maxit, tol);
semilogy(omegahist, 'b-');
set(ax16, 'YScale', 'log'); % this should be redundant
xlabel(ax16, "iteration"); ylabel(ax16, "\omega_k");
ylim(ax16,[1e-6,1.0])
legend(ax16,[string(thrs), "DP-BCG"] , 'Location', 'best');
grid(ax16, 'on');
hold(ax16, 'on');

figPath = fullfile(outDir, "dp_bf16.fig")
pdfPath = fullfile(outDir, "dp_bf16.pdf")
savefig(fig16, figPath);
exportgraphics(fig16, pdfPath, 'ContentType', 'vector');

%%%%%%%%%%%%%%%%%%%%%%%%
%m = 64
%%%%%%%%%%%%%%%%%%%%%%%%
m2 = 64;
X_ex2 = rand(n,m2);
b2 = A*X_ex2;
X_zero2 = zeros(n,m2);
fig64 = figure;
ax64 = axes(fig64);
hold(ax64,'on');

for i = 1:length(thrs)
	thr = thrs(i);
	linespec = linespecs(i);
	[~, ~, ~, omegahist] = bf_bcg(A, b2, applyM, X_ex2, maxit, tol, thr);
	semilogy(omegahist, linespec);
end

[~, omegahist] = PDPBCG(A, b2, X_zero2, applyM, X_ex2, maxit, tol);
semilogy(omegahist, 'b-');
set(ax64, 'YScale', 'log'); %this should be redundant
xlabel(ax64, "iteration"); ylabel(ax64, "\omega_k");
ylim(ax64,[1e-6,1.0]);
legend(ax64,[string(thrs), "DP-BCG"] , 'Location', 'best');
grid(ax64, 'on');
hold(ax64,'on');

figPath = fullfile(outDir, "dp_bf64.fig")
pdfPath = fullfile(outDir, "dp_bf64.pdf")
savefig(fig64,figPath);
exportgraphics(fig64, pdfPath, 'ContentType', 'vector');

