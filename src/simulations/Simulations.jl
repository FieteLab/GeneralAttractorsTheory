module Simulations
using LinearAlgebra: I, ⋅, mul!, norm  # import identity matrix and dot product
using Parameters: @with_kw_noshow
using Plots
using Term.Progress
using Statistics
using SparseArrays

import GeneralAttractors: save_simulation_history, moving_average, savepath
import GeneralAttractors: show_oneforms, show_oneforms!
import ..Can: AbstractCAN
import ..Manifolds: AbstractManifold, Manifoldℝ²

export Simulation, run_simulation
export Trajectory

include("trajectory.jl")
include("simulation.jl")
include("utils.jl")
end
