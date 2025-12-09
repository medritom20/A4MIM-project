function o_k = omega_error(A, x_k, x_ex, omegaInit)
% OMEGA_ERROR  Generaliized relative A-norm of the error.
%
% o_k = omega_error(A, x_k, x_ex, omegaDen)
%     A        SPD matrix.
%     x_k      current iterate.
%     x_ex     reference solution.
%     omegaDen optional precomputed trace(x_ex' * A * x_ex).

    if nargin < 4 || isempty(omegaInit)
        omegaInit = trace(x_ex' * A * x_ex);
        fprintf('Computed omegaInit = %.4e\n', omegaInit);
        fprintf('Profiling: Consider precomputing and passing as argument for efficiency\n');
    end

    err = x_ex - x_k;
    Aerr = A * err;
    num = sum(err .* Aerr, 'all'); % = trace(err' * A * err)
    o_k = sqrt(num / omegaInit);

    % Alternative (slower) implementation:
        % err = x_ex - x_k;
        % o_k = sqrt(trace(err' * A * err) / trace(x_ex' * A * x_ex));
end
