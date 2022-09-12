using GeneralAttractors

# specify simulation "chunks" were the input vector is constant
chunks = map(
    θ -> SimulationChunk([cos(θ), sin(θ)]; duration=600),
    [0, π/4, π/2, 3/4*π, π, -3/4*π, -π/2, -π/4]
) |> collect
chunks = SimulationChunk[SimulationChunk([0.0, 0.0]; duration=600), chunks...]

run_simulation(torus_attractor, chunks)
