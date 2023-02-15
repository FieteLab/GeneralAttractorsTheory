module GeneralAttractors
using LinearAlgebra: norm
using Plots
using Distances: PeriodicEuclidean, evaluate, UnionMetric, SphericalAngle, Euclidean
import Base.Iterators: product as ×  # cartesian product
using Term
using YAML, NPZ, JLD2
using ProjectSupervisor: Supervisor

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
include("animation.jl")

using .Simulations

import .Analysis


@info "Generating networks"
include("networks/mobius.jl")
include("networks/torus.jl")
include("networks/ring.jl")
include("networks/line.jl")
include("networks/sphere.jl")
include("networks/cylinder.jl")
include("networks/plane.jl")

export torus_maker, sphere_maker, ring_maker, mobius_maker, cylinder_maker, plane_maker, line_maker

end
