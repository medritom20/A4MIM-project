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
end
if nargin < 5 || isempty(Xtrue)
    error('PDPBCG requires Xtrue to evaluate omega.');
end
if nargin < 6 || isempty(maxIt)
    maxIt = 100;
    warning('PDPBCG:MaxItDefault', 'maxIt not provided; defaulting to %d.', maxIt);
end
if nargin < 7 || isempty(tol)
    tol = 1e-8;
    warning('PDPBCG:TolDefault', 'tol not provided; defaulting to %.1e.', tol);
end

omegaDen = sqrt(trace((Xtrue' * A) * Xtrue));
if omegaDen == 0
    error('The reference solution has zero A-norm; omega is undefined.');
end

X = X0;
R = B - A * X;
Z = applyM(R);
[P, ~] = qr(Z, 0);

omegaHist = zeros(maxIt + 1, 1);
E = Xtrue - X;
omega = sqrt(trace((E' * A) * E)) / omegaDen;
resTrace = trace(R' * R);
omegaHist(1) = omega;
fprintf('PDP-BCG: it.%3d\ttrace(r''*r) = %.2e\tomega = %.2e\n', 0, resTrace, omega);

for k = 1:maxIt
    AP = A * P;
    Gamma = (P' * AP) \ (P' * R);

    X = X + P * Gamma;
    R = R - AP * Gamma;

    Z = applyM(R);
    Delta = - (P' * AP) \ (P' * (A * Z));
    [P, ~] = qr(Z + P * Delta, 0);

    E = Xtrue - X;
    omega = sqrt(trace((E' * A) * E)) / omegaDen;
    resTrace = trace(R' * R);
    omegaHist(k + 1) = omega;
    fprintf('PDP-BCG: it.%3d\t\ttrace(r''*r) = %.2e\t\tomega = %.2e\n', k, resTrace, omega);

    if omega < tol
        omegaHist = omegaHist(1:k + 1);
        fprintf('PDP-BCG converged after %d iterations.\n', k);
        return
    end
end

omegaHist = omegaHist(1:maxIt + 1);
fprintf('PDP-BCG reached the maximum of %d iterations.\n', maxIt);

end
