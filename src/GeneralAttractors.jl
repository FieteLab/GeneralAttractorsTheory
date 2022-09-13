module GeneralAttractors
    include("io.jl")
    include("utils.jl")

    include("kernels.jl")
    include("can.jl")
    

    export CAN, Kernel
    export show_connectivity, plot_distance_function
    export ring_attractor, torus_attractor, mobius_attractor
    export MobiusEuclidean

    using .Kernels: AbstractKernel, Kernel, MexicanHatKernel, DiffOfExpKernel
    using .Can: CAN

    include("metrics.jl")
    include("viz.jl")
    include("networks.jl")

    include("simulations/Simulations.jl")

    using .Simulations
    export Simulation, ConstantChunk, run_simulation

end
 