module GeneralAttractors
using LinearAlgebra: norm
using Plots
using Measures
using Distances: PeriodicEuclidean, evaluate, UnionMetric, SphericalAngle
import Base.Iterators: product as ×  # cartesian product

include("io.jl")
include("utils.jl")

include("kernels.jl")
include("manifolds.jl")
include("can.jl")


export CAN, OneForm
export show_connectivity, plot_distance_function, show_oneforms, show_oneforms!

export load_simulation_history, save_data, load_data, save_model, load_model

using .Kernels
using .ManifoldUtils
using .Can: CAN, OneForm, offset_for_visual

include("metrics.jl")
include("viz.jl")
include("simulations/Simulations.jl")
include("analysis/Analysis.jl")

using .Simulations

import .Analysis

end
