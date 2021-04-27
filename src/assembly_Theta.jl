function assembly_Theta(X_domain, X_boundary, set_sigma)
    N_domain = size(X_domain)[1]
    N_boundary = size(X_boundary)[1]
    X_all = [X_domain;X_boundary]
    Theta = zeros((2*N_domain + N_boundary, 2*N_domain + N_boundary))
    
    for iter_i = 1:N_domain
        for iter_j = 1:N_domain
            Theta[iter_i,iter_j] = Delta_x_y_kappa(X_domain[iter_i,1], X_domain[iter_i,2], X_domain[iter_j,1], X_domain[iter_j,2], set_sigma)
        end
        # for loop is good!
        for iter_j = N_domain+1:2*N_domain+N_boundary
            Theta[iter_i,iter_j] = Delta_x_kappa(X_domain[iter_i,1], X_domain[iter_i,2], X_all[iter_j-N_domain,1], X_all[iter_j-N_domain,2], set_sigma)
            Theta[iter_j,iter_i] = Theta[iter_i,iter_j]
        end
    end
    for iter_i = N_domain+1:2*N_domain+N_boundary
        for iter_j = N_domain+1:2*N_domain+N_boundary
            Theta[iter_i,iter_j] = kappa(X_all[iter_i-N_domain,1], X_all[iter_i-N_domain,2], X_all[iter_j-N_domain,1], X_all[iter_j-N_domain,2], set_sigma)
        end
    end
    return Theta
end