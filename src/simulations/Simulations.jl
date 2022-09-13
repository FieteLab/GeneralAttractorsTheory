module Simulations
    using LinearAlgebra: I, ⋅, mul!, norm  # import identity matrix and dot product
    using Parameters: @with_kw_noshow
    using Plots
    using Term.Progress
    using Tullio
    using Statistics

    import GeneralAttractors: save_simulation_history, moving_average
    import ..Can: AbstractCAN

    export Simulation, ConstantChunk, run_simulation, RandomChunk

    include("chunk.jl")
    include("simulation.jl")
    include("utils.jl")
end