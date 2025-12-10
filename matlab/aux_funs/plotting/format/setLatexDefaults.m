function setLatexDefaults(fig)
%SETLATEXDEFAULTS  Use LaTeX interpreters in the given figure.
%   SETLATEXDEFAULTS()  uses current figure (gcf).

if nargin < 1 || isempty(fig)
    fig = gcf;
end

set(fig, ...
    'DefaultTextInterpreter',        'latex', ...
    'DefaultLegendInterpreter',      'latex', ...
    'DefaultAxesTickLabelInterpreter','latex');
end
