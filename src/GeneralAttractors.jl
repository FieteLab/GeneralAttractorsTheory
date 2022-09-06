module GeneralAttractors


    include("can.jl")

    export CAN, Kernel
    export show_connectivity

    using .Can: CAN, Kernel


    include("viz.jl")


end
 