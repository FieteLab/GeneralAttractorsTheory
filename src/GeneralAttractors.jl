module GeneralAttractors
using LinearAlgebra: norm
using Plots
using Distances: PeriodicEuclidean, evaluate, UnionMetric, SphericalAngle
import Base.Iterators: product as ×  # cartesian product

using BSON, YAML, NPZ

include("io.jl")
include("utils.jl")
include("embeddings.jl")

include("kernels.jl")
include("metrics.jl")
include("manifolds.jl")
include("can.jl")

export save_sim_data
export CAN, OneForm

using .Kernels
using .ManifoldUtils
using .Can: CAN, OneForm, offset_for_visual

include("viz.jl")
include("simulations/Simulations.jl")
include("analysis/Analysis.jl")

using .Simulations

import .Analysis


include("supervisor.jl")
using .ProjectSupervisor

end
