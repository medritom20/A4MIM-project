function [Xblk, omega_final, iter, omegahist] = bf_bcg(A, B, applyM, X_ex, maxit, tol, thr)
	%BF_BCG  BF-Omin / BF-ECG (Algorithm 2.3) 
	%
	%   [Xblk, omega_final, iter, omegahist] = bf_bcg(A, B, X_ex, maxit, tol, thr)
	%
	%   A      : SPD matrix (n x n)
	%   B      : the RHS (n x m)
	%   applyM : the preconditioner as a function handle
	%   X_ex   : the exact solution to the problem
	%   maxit  : maximum number of iterations
	%   tol    : stopping tolerance on \omega
	%   thr    : threshold of "almost linear dependence"
	%
	%   Outputs:
	%     Xblk         : n x m block of accumulated solution updates (X_k)
	%              The scalar solution is x = sum(Xblk, 2).
	%     omega_final  : final "A-residual"
	%     iter         : number of iterations performed
	%     omegahist    : history of \omega_k (length iter+1)
	% Solves the block system AX = B using the break free conjugate gradient method; no preconditioning, the initial point is zero.
	%--------------------------------------------------------------%


	% initial residual block and search directions are the block RHS
	R = B;                 % R_0 (n x t)
	Z = applyM(R);         % Z_1 = R_0

	[n, t] = size(R);

	% if the RHS is zero, the solution must be zero
	r0_vec = sum(R, 2);    % r_0
	normb  = norm(r0_vec);

	if normb == 0
		Xblk    = zeros(n, t)
		omega_final  = 0;
		iter    = 0;
		omegahist = 0;
		return;
	end

	% initialisation of the arrays
	omegahist       = zeros(maxit+1, 1);
	omegahist(1)    = 1.0;   
	Xblk          = zeros(n, t);
	fprintf('BF-BCG: it.%3d\trrank = %2d\tomega = %.2e\n', 0, size(Z,2), omegahist(1));


	fprintf('it  tcur  rrank   ||R||        ||Z||       omega     alpha_norm\n');
	for k = 1:maxit
		tcur = size(Z, 2);
		if tcur == 0
			warning('bf_ecg:allDropped', ...
				'All search directions dropped at iteration %d.', k);
			break;
		end

		% step 3: P_k = Z_k (Z_k^T A Z_k)^(-1/2)
		% C = Z_k^T A Z_K is a SPD matrix, so C^-1 is SPD and C^_1 = C^(-1/2) (C^(-1/2))^T, meaning we can compute the Cholesky factorization of C  = L L^T and C^(-1/2) = L^-T

		AZ = A * Z;                % n x tcur; we compute this only once
		C  = Z' * AZ;              % tcur x tcur

		% Cholesky: C = L * L^T
		L = chol(C, 'lower');      % assume A is SPD and basis is good enough

		% P_k and A P_k
		P  = Z / L';
		AP = AZ / L';

		% step 4: alpha_k = P_k^T R_{k-1}
		alpha = P' * R;            % tcur x t

		% steps 5–6:
		% X_k = X_{k-1} + P_k alpha_k
		% R_k = R_{k-1} - A P_k alpha_k
		Xblk = Xblk + P * alpha;
		R    = R    - AP * alpha;

		% update the residuum arrays, compute the norm of the combined residuum
		err = X_ex - Xblk; 
		omega = sqrt(trace(err' * A * err) / trace(X_ex' * A * X_ex));
		omegahist(k+1) = omega;

		if omega < tol
			break;
		end

		% steps 10–11: orthomin update
		%   beta_k = (A P_k)^T R_k
		%   Z_{k+1} = R_k - P_k beta_k
		beta = AP' * R;            % tcur x t
		Znew = R - P * beta;       % n x t
		Znew = applyM(Znew);

		% step 12: SVD based rank truncation
		%   (rank-revealing SVD done via MATLAB svd)
		[U, S, ~] = svd(Znew,0);
		svs = diag(S);
		fprintf('it=%3d  min SV = %.2e, max SV = %.2e, rrank=%d\n', k, min(svs), max(svs), rrank);
		relSV = svs/max(svs); % the paper computes relative singular values
		relSV = svs / max(svs);
		fprintf('rel min SV = %.2e, max rel SV = %.2e\n', min(relSV), max(relSV));
		rrank = sum(relSV > thr); %the rank of upper triangular matrix is the number of nonzero diagonal elements 
		Z = U(:,1:rrank); % Z_{k+1}; this assures Z has linearly independent columns

		% for diagnostics of stalling, we print plenty of information
		normR = norm(R, 'fro');
		normZ = norm(Znew, 'fro');
		alpha_norm = norm(alpha, 'fro');
		fprintf('%3d  %3d   %3d   %.2e   %.2e   %.2e   %.2e\n', ...
			k, tcur, rrank, normR, normZ, omega, alpha_norm);

		fprintf('%3d  %3d   %3d   %.2e   %.2e   %.2e   %.2e\n', ...
			k, tcur, rrank, normR, normZ, omega, alpha_norm);

		if rrank == 0
			warning('bf_bcg:allDroppedRRQR', ...
				'All directions dropped by RRQR at iteration %d.', k);
			break;
		end
	end

	iter    = k;
	omegahist = omegahist(1:iter+1);
	omega_final = omegahist(end);

end

