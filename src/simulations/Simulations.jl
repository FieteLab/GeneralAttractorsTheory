module Simulations
using LinearAlgebra: I, ⋅, mul!, norm  # import identity matrix and dot product
using Parameters: @with_kw_noshow
using Plots
using Term.Progress
using Statistics
using SparseArrays
import ForwardDiff: jacobian

import GeneralAttractors:
    save_simulation_history, moving_average, savepath, save_model, save_data, by_column
import GeneralAttractors: show_oneforms, show_oneforms!
import ..Can: AbstractCAN, offset_for_visual, OneForm, AbstractWeightOffset
import ..ManifoldUtils:
    AbstractManifold,
    Manifoldℝ²,
    Torus,
    Sphere,
    Mobius,
    sphere_embedding,
    ψx,
    ψy,
    ψz,
    ψ_t,
    ψ_θ1,
    ψ_θ2,
    CoverSpace,
    Ring

export Simulation, run_simulation
export Trajectory

include("trajectory.jl")
include("decoding.jl")
include("simulation.jl")
include("utils.jl")
end
