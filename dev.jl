using GeneralAttractors
using LinearAlgebra: I, ⋅  # import identity matrix and dot product
using BenchmarkTools: @benchmark
import GeneralAttractors.Can: AbstractCAN
using Plots

# TODO improve performance
# TODO check validity
# TODO add visualization

# TODO test with non-square A
# TODO test with 1D networks


# chunks = [
#     SimulationChunk([cos(0), sin(0)]),
#     SimulationChunk([cos(π/5), sin(π/5)]),
#     SimulationChunk([cos(π/2-π/5), sin(π/2-π/5)]),
# ]

chunks = map(
    θ -> SimulationChunk([cos(θ), sin(θ)]; duration=250),
    [0, π/4, π/2, 3/4*π, π, -3/4*π, -π/2, -π/4]
) |> collect


run_simulation(torus_attractor, chunks)
