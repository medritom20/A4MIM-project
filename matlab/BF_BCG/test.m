%%%%%%%%%%%%%%%%
% directory management
%%%%%%%%%%%%%%%%

srcDir = fullfile(pwd, "..");
bfDir = fullfile(srcDir, "BF_BCG");
pdpDir = fullfile(srcDir, "PDP_BCG");
dataDir = fullfile(srcDir, "data", "s3dkt3m2.mat");

% loading the data
A = load(dataDir).Problem.A;
n = length(A);

% iterations, tolerances
maxit = 200;
tol = 1e-10;
thrs = [1e-7, 1e-8];
linespecs = ["r:", "r-."];

%%%%%%%%%%%%%%%%
% preconditioner options
%%%%%%%%%%%%%%%%

opts.type = 'ict';        % threshold dropping
opts.droptol = 1e-5;      % drop tolerance
opts.diagcomp = 1e-2;     % global diagonal shift

L = ichol(A, opts);       
applyM = @(R) L' \ (L \ R);

%%%%%%%%%%%%%%%%
%m = 16
%%%%%%%%%%%%%%%%

m1 = 16;
X_ex1 = rand(n,m1);
b1 = A*X_ex1;
X_zero1 = zeros(n,m1);

fig = figure;
ax = axes(fig);
hold(ax,'on');

for i=1:length(thrs)
	thr = thrs(i);
	linespec = char(linespecs(i));
	[~, ~, ~, omegahist] = bf_bcg(A, b1, applyM, X_ex1, maxit, tol, thr);
	semilogy(ax,omegahist, linespec);
end

[~, omegahist] = PDPBCG(A, b1, X_zero1, applyM, X_ex1, maxit, tol);
semilogy(ax, omegahist, 'b-');
set(ax, 'YScale', 'log');
xlabel(ax,"iteration"); ylabel(ax,"\omega_k");
ylim(ax,[1e-6,1.0])
legend(ax,[string(thrs), "DP-BCG"], 'Location', 'best');
grid(ax,'on')
hold(ax,'off');

figPath = fullfile(outDir, "test.fig")
pdfPath = fullfile(outDir, "test.pdf")
savefig(fig,figPath);
exportgraphics(fig, pdfPath, 'ContentType', 'vector');


