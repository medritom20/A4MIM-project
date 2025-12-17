function [X, omegaHist] = PBFBCG(A, B, X0, applyM, Xtrue, maxIt, tol, svdThr)
% PBFBCG  Preconditioned BF-BCG with omega tracking.
%
% [X, omegaHist] = PBFBCG(A, B, X0, applyM, Xtrue, maxIt, tol)
%   A       ... SPD matrix
%   B       ... block right-hand side
%   X0      ... initial guess
%   applyM  ... @(u) M^{-1}u; pass [] for no preconditioner
%   Xtrue   ... reference solution (needed for omega)
%   maxIt   ... maximum number of iterations
%   tol     ... stopping tolerance on omega
%   svdThr  ... relative tolerance for singular-value filtering (default 1e-16)

if nargin < 4 || isempty(applyM)
    applyM = @(u) u;
    warning('PBFBCG:NoPreconditioner', 'No preconditioner provided; using identity.');
end
if nargin < 5 || isempty(Xtrue)
    error('PBFBCG requires Xtrue to evaluate omega.');
end
if nargin < 6 || isempty(maxIt)
    maxIt = 500;
    warning('PBFBCG:MaxItDefault', 'maxIt not provided; defaulting to %d.', maxIt);
end
if nargin < 7 || isempty(tol)
    tol = 1e-8;
    warning('PBFBCG:TolDefault', 'tol not provided; defaulting to %.1e.', tol);
end
if nargin < 8 || isempty(svdThr)
    svdThr = 1e-16;
    warning('PBFBCG:SvdTolDefault', 'svdThr not provided; defaulting to %.1e.', svdThr);
end

X = X0;
R = B - A * X;
Z = applyM(R);
[P, ~] = qr(Z, 0);

omegaInit = trace((Xtrue' * A) * Xtrue);                % precompute omegaInit
omegaHist = zeros(maxIt + 1, 1);
omegaHist(1) = omega_error(A, X, Xtrue, omegaInit);     % initial relative energy error; should be 1.0
fprintf('\n\tPDP-BCG: it.%3d\tomega = %.2e\n', 0, omegaHist(1));

for k = 1:maxIt
    AP = A * P;
    Gamma = (P' * AP) \ (P' * R);

    X = X + P * Gamma;
    R = R - AP * Gamma;

    Z = applyM(R);
    Delta = - (P' * AP) \ (P' * (A * Z));
    
    % DP-BCG: [P, ~] = qr(Z + P * Delta, 0);
    % BF-BCG: SVD with threshold
    [U, S, ~] = svd(Z + P * Delta, "econ");
        singVals = diag(S);
        sigmaMax = max(singVals);
        keep = singVals >= svdThr * sigmaMax;
    P = U(:, keep);

    omega = omega_error(A, X, Xtrue, omegaInit);
    omegaHist(k + 1) = omega;

        % Print iteration info
        fprintf('\tPBF-BCG: it.%3d\tomega = %.2e\t\tkappa_est(SVD) = %.2e', k, omega, sigmaMax / min(singVals(keep)));
        if any(~keep)
            fprintf('\t\t%2d/%2d singular values < %5.1f * (sigmaMax = %.3e) dropped; retaining %d.\n',  nnz(~keep), numel(singVals), svdThr, sigmaMax, nnz(keep));
        end
        fprintf('\n');
        
    if omega < tol
        omegaHist = omegaHist(1:k + 1);
        fprintf('PBF-BCG converged after %d iterations.\n', k);
        return
    end
end

omegaHist = omegaHist(1:maxIt + 1);
fprintf('PBF-BCG reached the maximum of %d iterations.\n', maxIt);

end
