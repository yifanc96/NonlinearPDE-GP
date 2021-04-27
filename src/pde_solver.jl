function PDE_solver(X_domain, X_boundary, max_iter, step_size, sol, L)
    N_domain = size(X_domain)[1]
    N_boundary = size(X_boundary)[1]
    rhs_f = zeros((N_domain,1))
    bdy_g = zeros((N_boundary,1))

    for iter_i = 1:N_domain
        rhs_f[iter_i] = f(X_domain[iter_i,1],X_domain[iter_i,2])
    end
    for iter_i = 1:N_boundary
        bdy_g[iter_i] = g(X_boundary[iter_i,1],X_boundary[iter_i,2])
    end
    
    J_array = zeros((max_iter+1,1))
    for iter_step = 1:max_iter
        temp_mtx = hcat(alpha*m*Diagonal((sol[:]).^(m-1)), Matrix(1.0I,N_domain,N_domain), zeros((N_domain,N_boundary)));
        temp_vec = [alpha*sol.^m-rhs_f; sol; bdy_g]
        J = L\temp_vec
        J = J'*J
        J_array[iter_step,1]=J[1]
        temp = L\(transpose(Matrix{Float64}(temp_mtx)))
        # J = 1/2*temp_vec'*Theta^{-1}*temp_vec
        grad_now = temp'*(L\temp_vec) #grad(temp_vec,sol)*Theta^{-1}*temp_vec
        GN_hessian_now = temp'*temp #grad(temp_vec,sol)*Theta^{-1}*grad(temp_vec,sol)
        sol = sol - step_size*GN_hessian_now\grad_now;
    end
    temp_vec=[alpha*sol.^m-rhs_f; sol; bdy_g]
    J = L\temp_vec
    J = J'*J
    J_array[max_iter+1,1]=J[1]
    return sol, J_array
end