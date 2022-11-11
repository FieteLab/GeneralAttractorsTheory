module Simulations
using LinearAlgebra: I, ⋅, mul!, norm  # import identity matrix and dot product
using Parameters: @with_kw_noshow
using Plots
using Term.Progress
using Statistics
using SparseArrays
import ForwardDiff: jacobian

import GeneralAttractors:
    save_simulation_history, moving_average, savepath, save_model, save_data
import GeneralAttractors: show_oneforms, show_oneforms!
import ..Can: AbstractCAN, offset_for_visual
import ..ManifoldUtils: AbstractManifold, Manifoldℝ², Torus, Sphere, sphere_embedding, ψx, ψy, ψz, ψxS², ψyS², ψzS²

export Simulation, run_simulation
export Trajectory

include("trajectory.jl")
include("simulation.jl")
include("utils.jl")
end
