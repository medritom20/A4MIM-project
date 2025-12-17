function defs = method_registry(context)
% METHOD_REGISTRY  Catalog of solver definitions.
%
% defs = METHOD_REGISTRY(context)
%   context can be either:
%       - results struct returned from driver (needs .precInfo)
%       - solver struct describing the experiment (needs .precCfg)
%
% The helper returns a column vector of style definitions with fields:
%   .name, .legend, .lineStyle, .color

if isfield(context, 'precInfo')
    precInfo = context.precInfo;
elseif isfield(context, 'precCfg')
    precInfo = context.precCfg;
elseif isfield(context, 'prec')
    precInfo = context.prec;
elseif isstruct(context) && isfield(context, 'type')
    precInfo = context;
else
    error('method_registry:missingPrec', ...
        'Context must provide precInfo/precCfg/prec.');
end

suffixPrec = legend_suffix(precInfo);

defs = [ ...
    struct('name', 'HS-BCG', 'legend', ['HS-BCG' suffixPrec], ...
           'lineStyle', ':',  'color', [0.8500 0.3250 0.0980], 'svdTol', []);
    struct('name', 'DP-BCG', 'legend', ['DP-BCG' suffixPrec], ...
           'lineStyle', '--', 'color', [0.0000 0.4470 0.7410], 'svdTol', []);
    struct('name', 'DR-BCG', 'legend', ['DR-BCG' suffixPrec], ...
           'lineStyle', '-',  'color', [0.0000 0.0000 0.0000], 'svdTol', []) ...
];

bfDefs = build_bfbcg_defs(context, suffixPrec);
defs = [defs; bfDefs];

defs = defs(:);
end

function bfDefs = build_bfbcg_defs(context, suffixPrec)
    bfDefs = struct('name', {}, 'legend', {}, 'lineStyle', {}, 'color', {}, 'svdTol', {});

    [svdTols, stylePref] = extract_bf_config(context);
    if isempty(svdTols)
        return;
    end

    lineStyles = resolve_bf_styles(stylePref, numel(svdTols));
    lineColors = resolve_bf_colors(context.colors, numel(svdTols) ) ;
    bfColor = [0.8500 0.1000 0.1000]; % force red hue

    bfDefs(numel(svdTols), 1) = struct('name', 'BF-BCG', 'legend', '', ...
                                       'lineStyle', '-.', 'color', bfColor, ...
                                       'svdTol', []);

    for k = 1:numel(svdTols)
        tol = svdTols(k);
        suffixTol = legend_svd_suffix(tol);
        bfDefs(k).name = 'BF-BCG'; 
        bfDefs(k).legend = sprintf('BF-BCG%s%s', suffixTol, suffixPrec);
        bfDefs(k).lineStyle = lineStyles{k};
        bfDefs(k).color = lineColors{k};
        bfDefs(k).svdTol = tol;
    end
end

function [svdTols, stylePref] = extract_bf_config(context)
    svdTols = [];
    stylePref = {};

    if isfield(context, 'svdTols') && ~isempty(context.svdTols)
        svdTols = context.svdTols(:)';
        if isfield(context, 'bfStyles') && ~isempty(context.bfStyles)
            styles = context.bfStyles;
            if iscell(styles)
                stylePref = styles(:)';
            end
        end
        return;
    end

    if isfield(context, 'sessions')
        buf = [];
        for s = 1:numel(context.sessions)
            runs = context.sessions(s).runs;
            for ir = 1:numel(runs)
                if strcmpi(runs(ir).method, 'BF-BCG') && ...
                        isfield(runs(ir).policy, 'svdTol') && ~isempty(runs(ir).policy.svdTol)
                    buf(end+1) = runs(ir).policy.svdTol; %#ok<AGROW>
                end
            end
        end
        if ~isempty(buf)
            svdTols = unique(buf, 'stable');
        end
    end
end

function styles = resolve_bf_styles(userStyles, n)
    palette = {':', '--', '-.', '-', ':'};
    styles = cell(1, n);
    for k = 1:n
        if k <= numel(userStyles) && ~isempty(userStyles{k})
            styles{k} = userStyles{k};
        else
            idx = min(k, numel(palette));
            styles{k} = palette{idx};
        end
    end
end

function linecolors = resolve_bf_colors(userStyles, n)
    defaultColor = 'red';
    linecolors = cell(1, n);
    for k = 1:n
        if k <= size(userStyles,1) && ~isempty(userStyles(k,:))
            linecolors{k} = userStyles(k,:);
        else
            linecolors{k} = defaultColor;
        end
    end
end

function suffix = legend_suffix(precInfo)
    if ~isfield(precInfo, 'type') || strcmpi(precInfo.type, 'none')
        suffix = '';
        return;
    end

    switch lower(precInfo.type)
        case 'ichol'
            [mant,exp] = format_scientific(precInfo.opts.droptol);
            if abs(mant - round(mant)) < 1e-10
                if abs(round(mant))-1 < 1e-10
                    tolTxt = sprintf('10^{%d}', exp);
                else
                    tolTxt = sprintf('%d \\times 10^{%d}', round(mant), exp);
                end
            else
                tolTxt = sprintf('%.1f \\times 10^{%d}', mant, exp);
            end
            suffix = sprintf(' + ICT($%s$)', tolTxt);
        otherwise
            suffix = sprintf(' + %s', upper(precInfo.type));
    end
end

function suffix = legend_svd_suffix(svdTol)
    [mant,exp] = format_scientific(svdTol);
    if abs(mant - round(mant)) < 1e-10
        if abs(round(mant))-1 < 1e-10
            tolTxt = sprintf('10^{%d}', exp);
        else
            tolTxt = sprintf('%d \\times 10^{%d}', round(mant), exp);
        end
    else
        tolTxt = sprintf('%.1f \\times 10^{%d}', mant, exp);
    end
    suffix = sprintf(' ($\\theta_{\\mathrm{SVD}} = %s$)', tolTxt);
end
