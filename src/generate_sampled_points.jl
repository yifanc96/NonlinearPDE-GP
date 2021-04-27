function sample_points(N_domain, N_boundary)
    # interior nodes
    X_domain = rand(Float64,(N_domain, 2))
    X_boundary = zeros((N_boundary, 2))
    
    # generate random boundary points on the faces of the domain 
    N_bd_each=convert(Int64, N_boundary/4)
    # bottom face
    X_boundary[1:N_bd_each, 1] = rand(Float64,(N_bd_each,1))
    # right face
    X_boundary[N_bd_each+1:2*N_bd_each, 1] .= 1
    X_boundary[N_bd_each+1:2*N_bd_each, 2] = rand(Float64,(N_bd_each,1))
    # top face
    X_boundary[2*N_bd_each+1:3*N_bd_each, 1] = rand(Float64,(N_bd_each,1))
    X_boundary[2*N_bd_each+1:3*N_bd_each, 2] .= 1
    # left face
    X_boundary[3*N_bd_each+1:N_boundary, 2] = rand(Float64,(N_bd_each,1))
    return X_domain, X_boundary
end