function o_k = omega_error(A,x_k,x_ex)

% Function used to calculate the special "A-norm" of the error. This norm
% is used in the numerical experiments section of the paper.
%
% A is positive definite matrix.
% x_k is the approximation of x_ex.

err = x_ex - x_k;
o_k = sqrt(trace(err' * A * err) / trace(x_ex' * A * x_ex));