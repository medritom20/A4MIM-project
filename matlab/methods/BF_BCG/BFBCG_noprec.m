function [Xblk, relres, iter, reshist] = BFBCG_noprec(A, B, tol, maxit)
%BF_BCG  BF-Omin / BF-ECG (Algorithm 2.3) 
%
%   [Xblk, relres, iter, reshist] = bf_ebg(A, B, tol, maxit)
%
%   A      : SPD matrix (n x n)
%   B      : the RHS (n x m)
%   tol    : stopping tolerance on ||sum(R_k)|| / ||sum(R_0)||
%   maxit  : maximum number of iterations
%
%   Outputs:
%     Xblk   : n x m block of accumulated solution updates (X_k)
%              The scalar solution is x = sum(Xblk, 2).
%     relres : final RELATIVE residual
%     iter   : number of iterations performed
%     reshist: history of ||r_k|| / ||r_0|| (length iter+1)
%--------------------------------------------------------------%


% Initial residual block and search directions are the block RHS
R = B;                 % R_0 (n x t)
Z = B;                 % Z_1 = R_0

[n, t] = size(R);

% If the RHS is zero, the solution must be zero
r0_vec = sum(R, 2);    % r_0
normb  = norm(r0_vec);

if normb == 0
    Xblk    = zeros(n, t);
    relres  = 0;
    iter    = 0;
    reshist = 0;
    return;
end

% Initialisation of the arrays
reshist       = zeros(maxit+1, 1);
reshist(1)    = 1.0;   % ||r_0|| / ||r_0|| = 1
Xblk          = zeros(n, t);

for k = 1:maxit

    tcur = size(Z, 2);
    if tcur == 0
        warning('bf_ecg:allDropped', ...
                'All search directions dropped at iteration %d.', k);
        break;
    end

    % ----------------------------------------
    % Step 3: P_k = Z_k (Z_k^T A Z_k)^(-1/2)
    % C = Z_k^T A Z_K is a SPD matrix, so C^-1 is SPD and C^_1 = C^(-1/2) (C^(-1/2))^T, meaning we can compute the Cholesky factorization of C  = L L^T and C^(-1/2) = L^-T

    % ----------------------------------------
    AZ = A * Z;                % n x tcur; we compute this only once
    C  = Z' * AZ;              % tcur x tcur

    % Cholesky: C = L * L^T
    L = chol(C, 'lower');      % assume A is SPD and basis is good enough

    % P_k and A P_k
    P  = Z / L';
    AP = AZ / L';

    % ----------------------------------------
    % Step 4: alpha_k = P_k^T R_{k-1}
    % ----------------------------------------
    alpha = P' * R;            % tcur x t

    % ----------------------------------------
    % Steps 5–6:
    % X_k = X_{k-1} + P_k alpha_k
    % R_k = R_{k-1} - A P_k alpha_k
    % ----------------------------------------
    Xblk = Xblk + P * alpha;
    R    = R    - AP * alpha;

    % Update the residuum arrays, compute the norm of the combined residuum
    r_vec  = sum(R, 2);
    relres = norm(r_vec) / normb;
    reshist(k+1) = relres;

    if relres < tol
        break;
    end

    % ----------------------------------------
    % Steps 10–11: Orthomin update
    %   beta_k = (A P_k)^T R_k
    %   Z_{k+1} = R_k - P_k beta_k
    % ----------------------------------------
    beta = AP' * R;            % tcur x t
    Znew = R - P * beta;       % n x t

    % ----------------------------------------
    % Step 12: RRQR(Z_{k+1}, sqrt(eps))
    %   (simple rank-revealing QR via MATLAB qr)
    % ----------------------------------------
    [Q, Rq] = qr(Znew, 0);     % Znew = Q * Rq
    diagR   = abs(diag(Rq)); %the magnitude of the diagonal elements of the upper triangular matrix Rq

    thr   = 1e-12 * max(diagR);    % this is our definition of "almost linearly dependent"
	    rrank = sum(diagR > thr); %the rank of upper triangular matrix is the number of nonzero diagonal elements 

    if rrank == 0
        warning('bf_ecg:allDroppedRRQR', ...
                'All directions dropped by RRQR at iteration %d.', k);
        break;
    end

    Z = Q(:, 1:rrank);         % Z_{k+1}; this assures Z has linearly independent columns
end

iter    = k;
reshist = reshist(1:iter+1);

end

