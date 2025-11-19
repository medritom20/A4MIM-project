function [omega_err, x_k] = DR_BCG_exp(A, b, x_0, x_ex, maxit)

x_k = x_0;
r_k = b - A * x_0;
[w_k,sig_k] = qr(r_k, 'econ');
s_k = w_k;

omega_err = zeros(maxit + 1, 1);
omega_err(1) = omega_error(A,x_0,x_ex);

for k = 1:maxit
    ksi_k_inv = s_k' * A * s_k;
    x_k = x_k + s_k * (ksi_k_inv \ sig_k);
    omega_err(k + 1) = omega_error(A, x_k, x_ex);
    [w_k,zeta] = qr(w_k - A * (s_k / ksi_k_inv), 'econ');
    s_k = w_k + s_k * zeta';
    sig_k = zeta * sig_k;
end