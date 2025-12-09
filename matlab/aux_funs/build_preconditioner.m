function [applyM, info] = build_preconditioner(A, precCfg)
% BUILD_PRECONDITIONER  Constructs a shared preconditioner description.
%
% [applyM, info] = build_preconditioner(A, precCfg)
%   applyM ... @(u) M^{-1} u
%   info   ... struct with fields: type, opts, L (if applicable), applyM

if nargin < 2 || isempty(precCfg) || ~isfield(precCfg, 'type')
    precCfg.type = 'none';
end

switch lower(precCfg.type)
    case 'ichol'
        opts = precCfg.opts;
        L = ichol(A, opts);
        applyM = @(x) L' \ (L \ x);
        info.type = 'ichol';
        info.opts = opts;
        info.L = L;

    case 'none'
        applyM = @(x) x;
        info.type = 'none';
        info.opts = struct();
        info.L = [];

    otherwise
        error('Unrecognized preconditioner type "%s".', precCfg.type);
end

info.applyM = applyM;

end
