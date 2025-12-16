function plot_convergence(omegaHist, solverLabel)

if nargin < 2 || isempty(solverLabel)
    solverLabel = 'PDP-BCG';
end

if isempty(omegaHist)
    warning('plot_convergence:noData', 'omega history is empty.');
    return;
end

k = 0:numel(omegaHist)-1;
hFig = figure('Name', 'PDP-BCG convergence');
ax = gca;
setLatexDefaults()
semilogy(k, omegaHist, 'LineWidth', 1.5);
grid on;

xlabel('$\mathrm{Iteration}$');
ylabel(['$\omega_k = \Bigl( ' ...
            '\frac{ ' ...
                '\textrm{tr} \left( ' ...
                    '(\textbf{x}-\textbf{x}_k)^{T} \textbf{A} (\textbf{x}-\textbf{x}_k)' ...
                '\right)' ...
            '}{' ...
                '\textrm{tr} \left( ' ...
                    '\textbf{x}^T \textbf{A} \textbf{x}' ...
                '\right)' ...
            '}' ...
        ' \Bigr)^{1/2}$']) ;
title('Convergence plot');

legend({solverLabel}, 'Location', 'best');
set(ax, 'YScale', 'log'); 

outDir = fullfile('.', 'out', 'results');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end
figPath = fullfile(outDir, 'omega_convergence.fig');
pdfPath = fullfile(outDir, 'omega_convergence.pdf');
savefig(hFig, figPath);
if verLessThan('matlab', '9.8')
    print(hFig, pdfPath, '-dpdf', '-vector');
else
    exportgraphics(hFig, pdfPath, 'ContentType', 'vector');
end

end
