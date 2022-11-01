module GeneralAttractors
    using LinearAlgebra: norm
    
    include("io.jl")
    include("utils.jl")

    include("kernels.jl")
    include("networks.jl")
    
    include("manifolds.jl")

    export CAN, Kernel
    export show_connectivity, plot_distance_function
    # export ring_attractor, torus_attractor, mobius_attractor, sphere_attractor
    export MobiusEuclidean
    export load_simulation_history, save_data, load_data, save_model, load_model

    using .Kernels
    using .Networks: CAN, IntegratorNetwork, VelocityNetwork

    include("metrics.jl")
    include("viz.jl")
    # include("simulations/Simulations.jl")
    # include("analysis/Analysis.jl")

    # using .Simulations
    # export Simulation, ConstantChunk, RandomChunk, run_simulation

    # import .Analysis

    # import .ManifoldUtils
    
    include("_networks.jl")

end
 