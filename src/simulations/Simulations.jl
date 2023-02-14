module Simulations
using LinearAlgebra: I, ⋅, mul!, norm  # import identity matrix and dot product
using Parameters: @with_kw_noshow
using Plots
using Term.Progress
using Statistics
using SparseArrays
import ForwardDiff: jacobian

import GeneralAttractors:
    moving_average,
    by_column,
    plot_can_vector_fields!
import GeneralAttractors: show_oneforms, show_oneforms!
import ..Can: AbstractCAN, offset_for_visual, OneForm, AbstractWeightOffset, SingleCAN
import ..ManifoldUtils:
    AbstractManifold,
    Manifoldℝ²,
    Torus,
    Sphere,
    Mobius,
    sphere_embedding,
    CoverSpace,
    Ring,
    apply_boundary_conditions!,
    AbstractVectorField

export Simulation, run_simulation
export Trajectory, ConstantTrajectory
export generate_groundtruth_data

include("trajectory.jl")
include("decoding.jl")
include("simulation.jl")
include("utils.jl")
include("single_can_groundtruth.jl")

end
