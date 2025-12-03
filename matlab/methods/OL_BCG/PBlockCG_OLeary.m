function [x_k1, res_norms] = PBlockCG_OLeary(A, b, X0, tol, maxiter, M_inv)
    r_0 = b - A * X0;              
    Z = apply_preconditioner(M_inv, r_0);  
    P_0 = Z;                       
    x_k1 = X0;
    res_norms = zeros(maxiter+1, 1);
    res_norms(1) = norm(r_0, 'fro');
    
    for k = 1:maxiter
        AP = A * P_0;
        gamma_k = (P_0' * AP) \ (Z' * r_0);     
        x_k1 = x_k1 + P_0 * gamma_k;                
        r_k = r_0 - AP * gamma_k;           
        res_norms(k+1) = norm(r_k, 'fro');
        
        if res_norms(k+1) < tol
            res_norms = res_norms(1:k+1);
            break;
        end
        
        Z_new = apply_preconditioner(M_inv, r_k);  
        Delta = (Z' * r_0) \ (Z_new' * r_k);          
        P_0 = Z_new + P_0 * Delta;                        
        
        
        r_0 = r_k;
        Z = Z_new;
    end
    
end

function Z = apply_preconditioner(M_inv, R)
    
    if isa(M_inv, 'function_handle')
        Z = zeros(size(R));
        for j = 1:size(R,2)
            Z(:,j) = M_inv(R(:,j));
        end
    else
        Z = M_inv * R;
    end
end