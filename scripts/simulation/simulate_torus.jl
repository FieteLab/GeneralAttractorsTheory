using GeneralAttractors
using Plots




simulation = Simulation(torus_attractor; η = 0.0)


MODE = :CONSTANT  # random or constant

if MODE == :CONSTANT
    # chunks = map(
    #     θ -> ConstantChunk([cos(θ), sin(θ)], simulation; duration=50),
    #     [0, π/4, π/2, 3/4*π, π, -3/4*π, -π/2, -π/4]
    # ) |> collect
    # chunks = ConstantChunk[ConstantChunk([0.0, 0.0], simulation; duration=100), chunks...]

    chunks = [
        ConstantChunk([0.0, 0.0], simulation; duration = 200),
        ConstantChunk([cos(π / 6), 0.0], simulation; duration = 200),
        ConstantChunk([cos(-π / 4), sin(-π / 4)], simulation; duration = 200),
        ConstantChunk([1.0, 0.0], simulation; duration = 800),
    ]
else
    chunks = [RandomChunk(simulation; duration = 150_000, μ₀ = 1.0, σ = 1)]
end
# plot(chunks[1])

h = @time run_simulation(
    simulation,
    chunks;
    frame_every_n = MODE == :CONSTANT ? 20 : nothing,
    discard_first_ms = MODE == :CONSTANT ? 0 : 5000,
    average_over_ms = 20,
)
