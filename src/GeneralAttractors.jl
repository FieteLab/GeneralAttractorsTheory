module GeneralAttractors
using LinearAlgebra: norm
using Plots
using Distances: PeriodicEuclidean, evaluate, UnionMetric, SphericalAngle
import Base.Iterators: product as ×  # cartesian product
using Term
using YAML, NPZ, JLD2

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
using .Can: CAN, OneForm, offset_for_visual, SingleCAN

include("viz.jl")
include("simulations/Simulations.jl")
include("analysis/Analysis.jl")

using .Simulations

import .Analysis


include("supervisor.jl")
using .ProjectSupervisor


include("networks/mobius.jl")
include("networks/torus.jl")
include("networks/ring.jl")
include("networks/sphere.jl")
export toruscan, toruscan_single, spherecan, ringcan, mobiuscan
export torus_maker


end
