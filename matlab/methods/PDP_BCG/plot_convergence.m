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
semilogy(k, omegaHist, 'LineWidth', 1.5);
grid on;

xlabel('$\mathrm{Iteration}$', 'Interpreter', 'latex');
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
        ' \Bigr)^{1/2}$'], ...
    'Interpreter', 'latex');
title('Convergence plot', 'Interpreter', 'latex');

set(gca, 'TickLabelInterpreter', 'latex');

legend({solverLabel}, 'Interpreter', 'latex', 'Location', 'best');

outDir = fullfile('.', 'out', 'results');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end
figPath = fullfile(outDir, 'omega_convergence.fig');
pdfPath = fullfile(outDir, 'omega_convergence.pdf');
savefig(hFig, figPath);
if verLessThan('matlab', '9.8')
    print(hFig, pdfPath, '-dpdf', '-painters');
else
    exportgraphics(hFig, pdfPath, 'ContentType', 'vector');
end

end
