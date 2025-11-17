function [omega_err, x] = Prec_DR_BCG_exp(A, b, x_0, x_ex, L, maxit)

x = x_0;
r = b - A * x_0;
[w,sig] = qr(L \ r);
s = L' \ w;

omega_err = zeros(maxit, 1);
omega_err(1) = omega_error(A,x_0,x_ex);
I = eye(size(s));

for k = 1:maxit
    ksi = (s' * A * s) \ I;
    x = x + s * ksi * sig;
    omega_err(k + 1) = omega_error(A, x, x_ex);
    [w,zeta] = qr(w - L \ A * s * ksi);
    s = L' \ w + s * zeta';
    sig = zeta * sig;
end