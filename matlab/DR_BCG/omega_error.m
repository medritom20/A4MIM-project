function o_k = omega_error(A,x_k,x_ex)

err = x_ex - x_k;
o_k = sqrt(trace((err)' * A * err) / trace(x_ex' * A * x_ex));