using GeneralAttractors
using Plots




simulation = Simulation(torus_attractor)


MODE = :RANDOM  # random or constant
if MODE == :CONSTANT
    chunks = map(
        θ -> ConstantChunk([cos(θ), sin(θ)], simulation; duration=250),
        [0, π/4, π/2, 3/4*π, π, -3/4*π, -π/2, -π/4]
    ) |> collect
    chunks = ConstantChunk[ConstantChunk([0.0, 0.0], simulation; duration=250), chunks...]
else
    chunks = [
        # ConstantChunk([0.0, 0.0], simulation; duration=500),
        RandomChunk(simulation; duration=100_000, μ₀=1.0, σ=1),
    ]
end
# plot(chunks[1])
 
h = run_simulation(simulation, chunks;  frame_every_n=nothing, discard_first_ms=1000, average_over_ms=15)
 

# TODO make it work
# TODO save gif in the right place