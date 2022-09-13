using GeneralAttractors

simulation = Simulation(torus_attractor)

# specify simulation "chunks" were the input vector is constant
chunks = map(
    θ -> ConstantChunk([cos(θ), sin(θ)], simulation; duration=250),
    [0, π/4, π/2, 3/4*π, π, -3/4*π, -π/2, -π/4]
) |> collect

# ! the first 0-v chunk is important to create the hex pattern
chunks = ConstantChunk[ConstantChunk([0.0, 0.0], simulation; duration=250), chunks...]

run_simulation(simulation, chunks)
 