module Simulations
    using LinearAlgebra: I, ⋅, mul!, norm  # import identity matrix and dot product
    using Parameters: @with_kw_noshow
    using Plots
    using Term.Progress
    using Tullio

    import ..Can: AbstractCAN

    export Simulation, SimulationChunk, run_simulation

    include("chunk.jl")
    include("simulation.jl")
    include("utils.jl")
end