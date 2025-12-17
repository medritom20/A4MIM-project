function [X, omegaHist] = PDPBCG(A, B, X0, applyM, Xtrue, maxIt, tol)
% PDPBCG  Preconditioned DP-BCG (Algorithm 8) with omega tracking.
%
% [X, omegaHist] = PDPBCG(A, B, X0, applyM, Xtrue, maxIt, tol)
%   A       ... SPD matrix
%   B       ... block right-hand side
%   X0      ... initial guess
%   applyM  ... @(u) M^{-1}u; pass [] for no preconditioner
%   Xtrue   ... reference solution (needed for omega)
%   maxIt   ... maximum number of iterations
%   tol     ... stopping tolerance on omega

if nargin < 4 || isempty(applyM)
    applyM = @(u) u;
    warning('PDPBCG:NoPreconditioner', 'No preconditioner provided; using identity.');
end
if nargin < 5 || isempty(Xtrue)
    error('PDPBCG requires Xtrue to evaluate omega.');
end
if nargin < 6 || isempty(maxIt)
    maxIt = 500;
    warning('PDPBCG:MaxItDefault', 'maxIt not provided; defaulting to %d.', maxIt);
end
if nargin < 7 || isempty(tol)
    tol = 1e-8;
    warning('PDPBCG:TolDefault', 'tol not provided; defaulting to %.1e.', tol);
end

X = X0;
R = B - A * X;
Z = applyM(R);
[P, ~] = qr(Z, 0);

omegaInit = trace((Xtrue' * A) * Xtrue);                % precompute omegaInit
omegaHist = zeros(maxIt + 1, 1);
omegaHist(1) = omega_error(A, X, Xtrue, omegaInit);     % initial relative energy error; should be 1.0
fprintf('\tPDP-BCG: it.%3d\t\tomega = %.2e\n', 0, omegaHist(1));

for k = 1:maxIt
    AP = A * P;
    Gamma = (P' * AP) \ (P' * R);

    X = X + P * Gamma;
    R = R - AP * Gamma;

    Z = applyM(R);
    Delta = - (P' * AP) \ (P' * (A * Z));
    [P, ~] = qr(Z + P * Delta, 0);

    omega = omega_error(A, X, Xtrue, omegaInit);
    omegaHist(k + 1) = omega;
    fprintf('\tPDP-BCG: it.%3d\t\tomega = %.2e\n', k, omega);
    if omega < tol
        omegaHist = omegaHist(1:k + 1);
        fprintf('PDP-BCG converged after %d iterations.\n', k);
        return
    end
end

omegaHist = omegaHist(1:maxIt + 1);
fprintf('PDP-BCG reached the maximum of %d iterations.\n', maxIt);

end
