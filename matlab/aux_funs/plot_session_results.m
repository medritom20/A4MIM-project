function plot_session_results(results, styleDefs, plotCfg)
% PLOT_SESSION_RESULTS  Visualizes omega histories from driver().
%
% plotCfg fields:
%   .layout       [rows cols]
%   .title        figure title (optional)
%   .outputDir    directory for saving (optional)
%   .figureSlug   base filename without extension (optional)
%   .saveOutputs  logical flag (default false)
%   .saveSingles  logical flag (default matches saveOutputs)
%   .experimentId experiment identifier (used for per-session figures)

if ~isfield(plotCfg, 'saveOutputs')
    plotCfg.saveOutputs = false;
end
if ~isfield(plotCfg, 'saveSingles')
    plotCfg.saveSingles = plotCfg.saveOutputs;
end
if ~isfield(plotCfg, 'experimentId')
    plotCfg.experimentId = [];
end

rows = plotCfg.layout(1);
cols = plotCfg.layout(2);

if plotCfg.saveOutputs && ~isempty(plotCfg.outputDir) && ~exist(plotCfg.outputDir, 'dir')
    mkdir(plotCfg.outputDir);
end

hFig = figure('Name', plotCfg.figureSlug, 'Color', 'w');
setLatexDefaults();
t = tiledlayout(rows, cols, 'TileSpacing', 'compact', 'Padding', 'compact');
if isfield(plotCfg, 'title') && ~isempty(plotCfg.title)
    title(t, plotCfg.title, 'Interpreter', 'tex');
end

sessions = results.sessions;
for isess = 1:numel(sessions)
    sess = sessions(isess);
    ax = nexttile(t);
    [kMax, legendEntries] = render_session(ax, sess, styleDefs);
    blockLabel = sprintf('$m = %d$', sess.blockSize);
    finalize_axes(ax, blockLabel, kMax, legendEntries, 'latex');

    if plotCfg.saveSingles && ~isempty(plotCfg.outputDir)
        singleFig = figure('Visible', 'off', 'Color', 'w');
        singleAx = axes('Parent', singleFig);
        [singleKMax, singleLegend] = render_session(singleAx, sess, styleDefs);
        singleTitle = format_single_title(plotCfg, sess);
        finalize_axes(singleAx, singleTitle, singleKMax, singleLegend, 'tex');
        singleSlug = build_single_slug(plotCfg, sess);
        save_figure(singleFig, plotCfg.outputDir, singleSlug);
        close(singleFig);
    end
end

if plotCfg.saveOutputs && ~isempty(plotCfg.outputDir)
    save_figure(hFig, plotCfg.outputDir, plotCfg.figureSlug);
end

end

function [kMax, legendEntries] = render_session(ax, session, styleDefs)
    hold(ax, 'on');
    grid(ax, 'on');
    set(ax, 'TickLabelInterpreter', 'latex');

    methods = session.runs;
    kMax = 0;
    legendEntries = cell(1, numel(methods));

    for im = 1:numel(methods)
        meth = methods(im);
        style = lookup_style(styleDefs, meth);
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
end

function finalize_axes(ax, titleStr, kMax, legendEntries, titleInterpreter)
    xlabel(ax, '$\mathrm{Iteration}$ $k$', 'Interpreter', 'latex');
    ylabel(ax, ['$\omega_k = \Bigl( ' ...
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
    title(ax, titleStr, 'Interpreter', titleInterpreter);

    set(ax, ...
        'YScale',               'log', ...
        'MinorGridLineStyle',   'none')

    if kMax > 0
        xlim(ax, [0, 100 * ceil(kMax / 100)]);
    end

    legend(ax, legendEntries, 'Interpreter', 'latex', 'Location', 'best');
end

function titleStr = format_single_title(plotCfg, session)
    parts = {};
    if ~isempty(plotCfg.experimentId)
        parts{end+1} = sprintf('Experiment %d', plotCfg.experimentId); %#ok<AGROW>
    end
    parts{end+1} = sprintf('rngSeed = %d', session.rngSeed); %#ok<AGROW>
    parts{end+1} = sprintf('m = %d', session.blockSize); %#ok<AGROW>
    titleStr = strjoin(parts, ' | ');
end

function slug = build_single_slug(plotCfg, session)
    tokens = {plotCfg.figureSlug};
    if ~isempty(plotCfg.experimentId)
        tokens{end+1} = sprintf('exp%d', plotCfg.experimentId); %#ok<AGROW>
    end
    tokens{end+1} = sprintf('seed%d', session.rngSeed); %#ok<AGROW>
    tokens{end+1} = sprintf('m%d', session.blockSize); %#ok<AGROW>
    slug = strjoin(tokens, '_');
end

function save_figure(figHandle, outputDir, slug)
    if isempty(outputDir) || isempty(slug)
        return;
    end
    figPath = fullfile(outputDir, [slug '.fig']);
    pdfPath = fullfile(outputDir, [slug '.pdf']);
    savefig(figHandle, figPath);
    exportgraphics(figHandle, pdfPath, 'ContentType', 'vector', 'BackgroundColor', 'white');
end

function style = lookup_style(styleDefs, runEntry)
methodName = '';
if isfield(runEntry, 'method') && ~isempty(runEntry.method)
    methodName = runEntry.method;
elseif isfield(runEntry, 'name') && ~isempty(runEntry.name)
    methodName = runEntry.name;
end
if isempty(methodName)
    methodName = 'unknown';
end
style = struct('color', [0.5 0.5 0.5], 'lineStyle', '-', 'legend', methodName);
runTol = [];
if isfield(runEntry, 'policy') && isfield(runEntry.policy, 'svdTol') && ~isempty(runEntry.policy.svdTol)
    runTol = runEntry.policy.svdTol;
end
for k = 1:numel(styleDefs)
    if strcmp(styleDefs(k).name, methodName)
        defTol = [];
        if isfield(styleDefs(k), 'svdTol')
            defTol = styleDefs(k).svdTol;
        end
        if ~isempty(defTol)
            if isempty(runTol)
                continue;
            end
            if abs(defTol - runTol) > max(1e-15, 1e-12 * abs(defTol))
                continue;
            end
        end
        style.color = styleDefs(k).color;
        style.lineStyle = styleDefs(k).lineStyle;
        style.legend = styleDefs(k).legend;
        return;
    end
end
end
