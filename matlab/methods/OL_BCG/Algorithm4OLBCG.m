function[x_k1] = Algorithm4OLBCG(A,b,x_0,tol,maxiter)
    r_0 = b - A * x_0;
    [~,numCols] = size(b);
    Phi_0 = eye(numCols);
    P_0 = r_0 * Phi_0;
    iter = 0;
    while norm(r_0,'fro') > tol && iter < maxiter

        gamma_k = (P_0.'*A*P_0) \ ((r_0.') * r_0);
        x_k1 = x_0 + P_0 * gamma_k;
        r_k = r_0 - A * (P_0 * gamma_k);
        delta_k = ((r_0.') * r_0) \ ((r_k.') * r_k);
        P_k = (r_k + P_0 * delta_k); 
        

        x_0 = x_k1;
        r_0 = r_k;
        P_0 = P_k;
        iter = iter + 1;
    end
end