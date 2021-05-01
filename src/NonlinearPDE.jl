module NonlinearPDE

# Write your package code here.
    import LinearAlgebra: Diagonal, I
    include("truth_solution.jl")
    export u, f, g
    include("kernel_functions.jl")
    export kappa, Delta_x_kappa, Delta_x_y_kappa
    include("generate_sampled_points.jl")
    export sample_points
    include("assembly_Theta.jl")
    export assembly_Theta
    include("pde_solver.jl")
    export PDE_solver
    include("extend_solution.jl")
    export construct_inpTheta, extend_solution, get_extended_solution
    include("error_calculation.jl")
    export solution_error
end
