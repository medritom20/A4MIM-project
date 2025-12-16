function [X, flag, relres, iter, resHist, omegaHist] = HS_PBCG(A, B, Mfun, X0, tol, maxit, phiType, Xtrue)
% HS_PBCG  Preconditioned HS-BCG (Algorithm 4 with Phi_k control)
%
% [X, flag, relres, iter, resHist, omegaHist] = HS_PBCG(A, B, Mfun, X0, tol, maxit, phiType, Xtrue)
%     A        SPD matrix.
%     B        right-hand side block (n x m).
%     Mfun     preconditioner handle or matrix that applies M^{-1}.
%     X0       initial guess.
%     tol      stopping tolerance on ||R_k||_F / ||R_0||_F.
%     maxit    maximum number of iterations.
%     phiType  scaling applied to search directions ('eye', 'chol', 'a-orth', or function handle).
%     Xtrue    optional reference solution for omega history.
%
% Outputs mirror MATLAB conventions: FLAG = 0 converged, 1 hit maxit, 2 breakdown.

if nargin < 7 || isempty(phiType)
    phiType = 'eye';
end
if nargin < 8
    Xtrue = [];
end

n = size(A, 1);
if size(B, 1) ~= n
    error('HS_PBCG:dimensionMismatch', 'A and B dimensions are inconsistent.');
end

X = X0;
R = B - A * X;
Z = apply_preconditioner(Mfun, R);
P = scale_block(A, Z, phiType, 0);

normR0 = norm(R, 'fro');
flag = 1;
iter = 0;
if normR0 == 0
    relres = 0;
    resHist = 0;
    if ~isempty(Xtrue)
        omegaHist = omega_error(A, X, Xtrue);
    else
        omegaHist = [];
    end
    flag = 0;
    return;
end

resHist = zeros(maxit + 1, 1);
resHist(1) = 1;
trackOmega = nargin >= 8 && ~isempty(Xtrue);
if trackOmega
    omegaDen = trace((Xtrue' * A) * Xtrue);
    omegaHist = zeros(maxit + 1, 1);
    omegaHist(1) = omega_error(A, X, Xtrue, omegaDen);
else
    omegaHist = [];
    omegaDen = [];
end
omegaPrint = NaN;
if trackOmega
    omegaPrint = omegaHist(1);
end
resTrace = trace(R' * R);
fprintf('HS-BCG: it.%3d\ttrace(r''*r) = %.2e\tomega = %.2e\n', 0, resTrace, omegaPrint);

for k = 1:maxit
    AP = A * P;
    GramP = P' * AP;
    RHS = Z' * R;

    [Gamma, gammaBroken] = solve_sym_pos(GramP, RHS, 'gamma', k);
    if gammaBroken
        flag = 2;
        iter = k - 1;
        break;
    end

    X = X + P * Gamma;
    Rnew = R - AP * Gamma;
    resNorm = norm(Rnew, 'fro');
    relres = resNorm / normR0;
    resHist(k + 1) = relres;
    iter = k;
    if trackOmega
        omegaHist(k + 1) = omega_error(A, X, Xtrue, omegaDen);
        omegaPrint = omegaHist(k + 1);
    else
        omegaPrint = NaN;
    end
    resTrace = trace(Rnew' * Rnew);
    fprintf('HS-BCG: it.%3d\ttrace(r''*r) = %.2e\tomega = %.2e\n', k, resTrace, omegaPrint);

    if relres <= tol
        flag = 0;
        break;
    end

    Znew = apply_preconditioner(Mfun, Rnew);
    [Delta, deltaBroken] = solve_linear(RHS, Znew' * Rnew, 'delta', k);
    if deltaBroken
        flag = 2;
        % break;
    end

    Praw = Znew + P * Delta;
    P = scale_block(A, Praw, phiType, k);

    R = Rnew;
    Z = Znew;
end

relres = resHist(iter + 1);
resHist = resHist(1:iter + 1);
if trackOmega
    omegaHist = omegaHist(1:iter + 1);
end

end

function Z = apply_preconditioner(Mfun, R)
if isa(Mfun, 'function_handle')
    Z = zeros(size(R));
    for j = 1:size(R, 2)
        Z(:, j) = Mfun(R(:, j));
    end
else
    Z = Mfun * R;
end
end

function P = scale_block(A, block, phiType, iter)
if isa(phiType, 'function_handle')
    Phi = phiType(A, block, iter);
    P = block * Phi;
    return;
end

switch lower(phiType)
    case 'eye'
        P = block;
    case 'chol'
        P = block / chol_safe(block' * block, 'chol', iter);
    case {'a-orth', 'chol-a'}
        Gram = block' * (A * block);
        P = block / chol_safe(Gram, 'chol-a', iter);
    otherwise
        error('HS_PBCG:unknownPhi', 'Unsupported phiType "%s".', phiType);
end
end

function U = chol_safe(G, tag, iter)
G = (G + G') * 0.5;
[U, p] = chol(G);
if p ~= 0
    warning('HS_PBCG:phiFallback', ...
        'Cholesky failed for phiType "%s" at iteration %d; reverting to identity.', tag, iter);
    U = eye(size(G, 1));
end
end

function [Xsol, broken] = solve_sym_pos(M, RHS, tag, iter)
broken = false;
M = (M + M') * 0.5;
[L, p] = chol(M, 'lower');
if p ~= 0
    warning('HS_PBCG:GammaBreakdown', ...
        'Matrix for %s at iteration %d is singular.', tag, iter);
    Xsol = zeros(size(M, 1), size(RHS, 2));
    broken = true;
    return;
end
Y = L \ RHS;
Xsol = L' \ Y;
end

function [Xsol, broken] = solve_linear(A, B, tag, iter)
broken = false;
if rcond(A) < 1e-14
    warning('HS_PBCG:DeltaBreakdown', ...
        'Coefficient matrix for %s near singular at iteration %d.', tag, iter);
    Xsol = zeros(size(A, 2), size(B, 2));
    broken = true;
    return;
end
Xsol = A \ B;
end
