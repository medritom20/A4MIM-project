function [X, omegaHist] = BFBCG(A, B, X0, applyM, Xtrue, maxIt, tol, svdTol)
% BFBCG  Breakdown-free BCG variant based on PDP-BCG with SVD filtering.
%
% [X, omegaHist] = BFBCG(A, B, X0, applyM, Xtrue, maxIt, tol, svdTol)
%   A        ... SPD matrix
%   B        ... block right-hand side
%   X0       ... initial guess
%   applyM   ... @(u) M^{-1}u; pass [] for no preconditioner
%   Xtrue    ... reference solution (needed for omega)
%   maxIt    ... maximum number of iterations
%   tol      ... stopping tolerance on omega
%   svdTol   ... relative tolerance for singular-value filtering (default 1e-10)

if nargin < 4 || isempty(applyM)
    applyM = @(u) u;
end
if ~isa(applyM, 'function_handle')
    M = applyM;
    applyM = @(x) M * x;
end
if nargin < 5 || isempty(Xtrue)
    error('BFBCG requires Xtrue to evaluate omega.');
end
if nargin < 6 || isempty(maxIt)
    maxIt = 500;
    warning('BFBCG:MaxItDefault', 'maxIt not provided; defaulting to %d.', maxIt);
end
if nargin < 7 || isempty(tol)
    tol = 1e-8;
    warning('BFBCG:TolDefault', 'tol not provided; defaulting to %.1e.', tol);
end
if nargin < 8 || isempty(svdTol)
    svdTol = 1e-10;
end

X = X0;
R = B - A * X;
Z = apply_block_preconditioner(applyM, R);
[P, ~] = qr(Z, 0);

omegaInit = trace((Xtrue' * A) * Xtrue);
omegaHist = zeros(maxIt + 1, 1);
omegaHist(1) = omega_error(A, X, Xtrue, omegaInit);
fprintf('\tBF-BCG: it.%3d\tomega = %.2e\n', 0, omegaHist(1));

for k = 1:maxIt
    if isempty(P)
        warning('BFBCG:emptyBasis', 'P became empty at iteration %d.', k);
        omegaHist = omegaHist(1:k);
        return;
    end

    AP = A * P;
    C = P' * AP;
    [pinvC, dropped] = filtered_pinv(C, svdTol);

    if dropped > 0
        fprintf('\tBF-BCG: dropped %d singular directions (rank %d/%d) at it.%d\n', ...
                dropped, size(C, 1) - dropped, size(C, 1), k - 1);
    end

    Gamma = pinvC * (P' * R);

    X = X + P * Gamma;
    R = R - AP * Gamma;

    Z = apply_block_preconditioner(applyM, R);
    AZ = A * Z;
    Delta = -pinvC * (P' * AZ);
    [P, ~] = qr(Z + P * Delta, 0);

    omega = omega_error(A, X, Xtrue, omegaInit);
    omegaHist(k + 1) = omega;
    fprintf('\tBF-BCG: it.%3d\tomega = %.2e\n', k, omega);
    if omega < tol
        omegaHist = omegaHist(1:k + 1);
        fprintf('BF-BCG converged after %d iterations.\n', k);
        return
    end
end

omegaHist = omegaHist(1:maxIt + 1);
fprintf('BF-BCG reached the maximum of %d iterations.\n', maxIt);

end

function [pinvC, dropped] = filtered_pinv(C, svdTol)
    if isempty(C)
        pinvC = zeros(size(C'));
        dropped = 0;
        return;
    end

    [U, S, V] = svd(C, 'econ');
    singVals = diag(S);

    if isempty(singVals)
        pinvC = zeros(size(C'));
        dropped = 0;
        return;
    end

    sigmaMax = max(singVals);
    if sigmaMax == 0
        keep = false(size(singVals));
        keep(1) = true;
    else
        keep = singVals >= svdTol * sigmaMax;
        if ~any(keep)
            [~, idx] = max(singVals);
            keep(idx) = true;
        end
    end

    Sinv = diag(1 ./ singVals(keep));
    pinvC = V(:, keep) * (Sinv * U(:, keep)');
    dropped = numel(singVals) - nnz(keep);
end

function Z = apply_block_preconditioner(op, R)
    if isa(op, 'function_handle')
        try
            Z = op(R);
        catch
            Z = zeros(size(R), 'like', R);
            for j = 1:size(R, 2)
                Z(:, j) = op(R(:, j));
            end
        end
    else
        Z = op * R;
    end
end
