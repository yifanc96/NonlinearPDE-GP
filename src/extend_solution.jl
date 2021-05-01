function construct_inpTheta(X_test, X_domain, X_boundary, set_sigma)
    N_test = size(X_test)[1]
    N_domain = size(X_domain)[1]
    N_boundary = size(X_boundary)[1]
    inpTheta = zeros((N_test, 2*N_domain+N_boundary))
    for iter_i = 1:N_test
        for iter_j = 1:N_domain
            inpTheta[iter_i,iter_j] = Delta_x_kappa(X_test[iter_i,1],X_test[iter_i,2],X_domain[iter_j,1],X_domain[iter_j,2],set_sigma)
        end
        for iter_j = N_domain+1:2*N_domain
            inpTheta[iter_i,iter_j] = kappa(X_test[iter_i,1],X_test[iter_i,2],X_domain[iter_j-N_domain,1],X_domain[iter_j-N_domain,2],set_sigma)
        end
        for iter_j = 2*N_domain+1:2*N_domain+N_boundary
            inpTheta[iter_i,iter_j] = kappa(X_test[iter_i,1],X_test[iter_i,2],X_boundary[iter_j-2*N_domain,1],X_boundary[iter_j-2*N_domain,2],set_sigma)
        end
    end
    return inpTheta
end

function extend_solution(X_test, X_domain, X_boundary, L, v, set_sigma)
    inp_Theta = construct_inpTheta(X_test, X_domain, X_boundary, set_sigma)
    return inp_Theta*(L'\(L\v))
end

function get_extended_solution(num_pts, X_domain, X_boundary, L, sol, set_sigma)
    x = range(0,1,length=num_pts)
    y = range(0,1,length=num_pts)
    xx = x' .* ones(num_pts)
    yy = ones(num_pts)' .* y
    
    xxv = Matrix{Float64}(reshape(xx,(num_pts^2,1)))
    yyv = Matrix{Float64}(reshape(yy,(num_pts^2,1)))

    X_test = hcat(xxv,yyv)

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

#     X_test=hcat([[x[i];y[j]] for i = 1 : num_pts for j = 1 : num_pts])
    temp_vec = [alpha*sol.^m-rhs_f; sol; bdy_g]
    u_extended = Matrix{Float64}(reshape(extend_solution(X_test, X_domain, X_boundary, L, temp_vec, set_sigma),(num_pts,num_pts)))
    return x, y, u_extended
end