module GeneralAttractors
using LinearAlgebra: norm
using Plots
using Distances: PeriodicEuclidean, evaluate, UnionMetric, SphericalAngle, Euclidean
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
include("animation.jl")

using .Simulations

import .Analysis


include("supervisor.jl")
using .ProjectSupervisor

@info "Generating networks"
include("networks/mobius.jl")
include("networks/torus.jl")
include("networks/ring.jl")
include("networks/sphere.jl")
include("networks/cylinder.jl")
include("networks/plane.jl")

export toruscan, toruscan_single, 
        spherecan, spherecan_single, 
        ringcan, ringcan_single, 
        mobiuscan, mobiuscan_single, 
        cylindercan, cylindercan_single,
        planecan, planecan_single
export torus_maker, sphere_maker, ring_maker, mobius_maker, cylinder_maker, plane_maker

end
