module GeneralAttractors


    include("can.jl")

    export CAN, Kernel
    export show_connectivity, plot_distance_function
    export ring_attractor, torus_attractor
    export MobiusEuclidean

    using .Can: CAN, Kernel

    include("metrics.jl")
    include("viz.jl")
    include("networks.jl")


end
 