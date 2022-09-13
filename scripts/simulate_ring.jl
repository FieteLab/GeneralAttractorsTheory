using GeneralAttractors

simulation = Simulation(
    ring_attractor; frame_every_n=10, η=0.1, b₀=1.0
)

chunks = map(
    θ -> SimulationChunk([Float64(θ)], simulation; duration=250),
    [0, -1, 0, 1]
) |> collect

run_simulation(simulation, chunks)
