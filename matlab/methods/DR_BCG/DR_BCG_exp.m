function [omega_err, x_k] = DR_BCG_exp(A, b, x_0, x_ex, maxit)

% Function meant for the experiments. Clculates only the special "A-norm"
% of the error. Also returns the last iterations x_k.
%
% A is a positive definite matrix.
% b is a matric right hand side.
% x_ex is an exact matrix solution to a problem A * x = b.
% x_0 is initial guess.
% maxit is the number of iterations.

    x_k = x_0;
    r_k = b - A * x_0;
    [w_k,sig_k] = qr(r_k, 'econ');
    s_k = w_k;

    omegaInit = trace((x_ex' * A) * x_ex);                  % precompute omegaInit
    omega_err = zeros(maxit + 1, 1);
    omega_err(1) = omega_error(A, x_0, x_ex, omegaInit);    % initial relative error; should be 1.0
        fprintf('\tDR-BCG:\tit.%3d\tomega = %.2e\n', 0, omega_err(1));
    for k = 1:maxit
        ksi_k_inv = s_k' * A * s_k;
        x_k = x_k + s_k * (ksi_k_inv \ sig_k);
        omega_err(k + 1) = omega_error(A, x_k, x_ex, omegaInit);
            fprintf('\tDR-BCG:\tit.%3d\tomega = %.2e\n', k, omega_err(k + 1));
        [w_k,zeta] = qr(w_k - A * (s_k / ksi_k_inv), 'econ');
        s_k = w_k + s_k * zeta';
        sig_k = zeta * sig_k;
    end

end