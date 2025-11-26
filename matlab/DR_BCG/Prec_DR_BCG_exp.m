function [omega_err, x_k] = Prec_DR_BCG_exp(A, b, x_0, x_ex, L, maxit)

% Function meant for the experiments. Clculates only the special "A-norm"
% of the error. Also returns the last iterations x_k.
%
% A is a positive definite matrix.
% b is a matric right hand side.
% x_ex is an exact matrix solution to a problem A * x = b.
% x_0 is initial guess.
% L is the cholesky factoritzations of the preconditioner M = L * L'
% maxit is the number of iterations.

x_k = x_0;
r_k = b - A * x_0;
[w_k,sig_k] = qr(L \ r_k, 'econ');
s_k = L' \ w_k;

omega_err = zeros(maxit + 1, 1);
omega_err(1) = omega_error(A,x_0,x_ex);

for k = 1:maxit
    ksi_k_inv = s_k' * A * s_k;
    x_k = x_k + s_k * (ksi_k_inv \ sig_k);
    omega_err(k + 1) = omega_error(A, x_k, x_ex);
    [w_k,zeta] = qr(w_k - L \ (A * (s_k / ksi_k_inv)), 'econ');
    s_k = L' \ w_k + s_k * zeta';
    sig_k = zeta * sig_k;
end