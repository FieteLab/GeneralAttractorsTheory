using GeneralAttractors

simulation = Simulation(
    ring_attractor; η=0.1, b₀=1.0
)


chunks = map(
    θ -> ConstantChunk([Float64(θ)], simulation; duration=250),
    [0, -1, 0, 1]
) |> collect

run_simulation(simulation, chunks; frame_every_n=10)
# TODO fix plotting for this