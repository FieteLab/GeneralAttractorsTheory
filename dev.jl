using GeneralAttractors
import GeneralAttractors.Can: CAN
using GeneralAttractors.Kernels
using Distances

# TODO test with non-square A
# TODO test with 1D networks


simulation = Simulation(
    ring_attractor
)


chunks = map(
    θ -> SimulationChunk([cos(θ)]; duration=600),
    [0, π/4, π/2, 3/4*π, π, -3/4*π, -π/2, -π/4]
) |> collect
chunks = SimulationChunk[SimulationChunk([10.0, 0.0], simulation; duration=250)]

run_simulation(can, chunks)
