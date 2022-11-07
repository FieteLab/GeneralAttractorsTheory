module GeneralAttractors
using LinearAlgebra: norm

include("io.jl")
include("utils.jl")

include("kernels.jl")
include("manifolds.jl")
include("can.jl")


export CAN, Kernel
export show_connectivity, plot_distance_function,  show_oneforms
# export ring_attractor, torus_attractor, mobius_attractor, sphere_attractor
export torus_attractor

export load_simulation_history, save_data, load_data, save_model, load_model

using .Kernels
using .Manifolds
using .Can: CAN, OneForm

include("metrics.jl")
include("viz.jl")
include("simulations/Simulations.jl")
include("analysis/Analysis.jl")

using .Simulations

import .Analysis

include("networks.jl")

end
