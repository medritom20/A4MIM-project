function plot_session_results(results, styleDefs, plotCfg)
% PLOT_SESSION_RESULTS  Visualizes omega histories from driver().
%
% plotCfg fields:
%   .layout      [rows cols]
%   .title       figure title (optional)
%   .outputDir   directory for saving (optional)
%   .figureSlug  base filename without extension (optional)
%   .saveOutputs logical flag (default false)

% arguments
%     results struct
%     styleDefs struct
%     plotCfg.layout (1,2) double
%     plotCfg.title char = ''
%     plotCfg.outputDir char = ''
%     plotCfg.figureSlug char = 'driver_plot'
%     plotCfg.saveOutputs logical = false
% end

rows = plotCfg.layout(1);
cols = plotCfg.layout(2);

hFig = figure('Name', plotCfg.figureSlug, 'Color', 'w');
setLatexDefaults();
t = tiledlayout(rows, cols, 'TileSpacing', 'compact', 'Padding', 'compact');
if ~isempty(plotCfg.title)
    title(t, plotCfg.title, 'Interpreter', 'tex');
end

sessions = results.sessions;
for isess = 1:numel(sessions)
    ax = nexttile(t);
    hold(ax, 'on');
    grid(ax, 'on');
    set(ax, 'TickLabelInterpreter', 'latex');

    blockLabel = sprintf('$m = %d$', sessions(isess).blockSize);
    methods = sessions(isess).runs;
    kMax = 0;
    legendEntries = cell(1, numel(methods));

    for im = 1:numel(methods)
        meth = methods(im);
        style = lookup_style(styleDefs, meth.name);
        if isfield(meth, 'style') && ~isempty(fieldnames(meth.style))
            if isfield(meth.style, 'color'),     style.color = meth.style.color; end
            if isfield(meth.style, 'lineStyle'), style.lineStyle = meth.style.lineStyle; end
            if isfield(meth.style, 'lineWidth'), style.lineWidth = meth.style.lineWidth; end
        end
        if ~isfield(style, 'lineWidth') || isempty(style.lineWidth)
            style.lineWidth = 1.5;
        end

        omegaHist = meth.results.omega;
        k = 0:numel(omegaHist) - 1;
        if ~isempty(k)
            kMax = max(kMax, k(end));
        end
        semilogy(ax, k, omegaHist, ...
            'Color', style.color, ...
            'LineStyle', style.lineStyle, ...
            'LineWidth', style.lineWidth);

        if isfield(meth, 'label') && ~isempty(meth.label)
            legendEntries{im} = meth.label; %#ok<AGROW>
        elseif isfield(style, 'legend')
            legendEntries{im} = style.legend; %#ok<AGROW>
        else
            legendEntries{im} = meth.name; %#ok<AGROW>
        end
    end

    xlabel('$\mathrm{Iteration}$ $k$');
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
            ' \Bigr)^{1/2}$'] );
    title(ax, blockLabel);

    set(ax, ...
        'YScale',               'log', ...
        'MinorGridLineStyle',   'none')

    if kMax > 0
        xlim(ax, [0, 100 * ceil(kMax / 100)]);
    end

    legend(ax, legendEntries, 'Interpreter', 'latex', 'Location', 'best');
end

if plotCfg.saveOutputs && ~isempty(plotCfg.outputDir)
    if ~exist(plotCfg.outputDir, 'dir')
        mkdir(plotCfg.outputDir);
    end
    figPath = fullfile(plotCfg.outputDir, [plotCfg.figureSlug '.fig']);
    pdfPath = fullfile(plotCfg.outputDir, [plotCfg.figureSlug '.pdf']);
    savefig(hFig, figPath);
    exportgraphics(hFig, pdfPath, 'ContentType', 'vector', 'BackgroundColor', 'white');
end

end

function style = lookup_style(styleDefs, methodName)
style = struct('color', [0.5 0.5 0.5], 'lineStyle', '-', 'legend', methodName);
for k = 1:numel(styleDefs)
    if strcmp(styleDefs(k).name, methodName)
        style.color = styleDefs(k).color;
        style.lineStyle = styleDefs(k).lineStyle;
        style.legend = styleDefs(k).legend;
        return;
    end
end
end
