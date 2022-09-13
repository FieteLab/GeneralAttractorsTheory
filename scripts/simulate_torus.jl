using GeneralAttractors

simulation = Simulation(torus_attractor)

# specify simulation "chunks" were the input vector is constant


MODE = :RANDOM  # random or constant

if MODE == :CONSTANT
    chunks = map(
        θ -> ConstantChunk([cos(θ), sin(θ)], simulation; duration=250),
        [0, π/4, π/2, 3/4*π, π, -3/4*π, -π/2, -π/4]
    ) |> collect
    chunks = ConstantChunk[ConstantChunk([0.0, 0.0], simulation; duration=250), chunks...]

else
    chunks = [
        ConstantChunk([0.0, 0.0], simulation; duration=500),
        RandomChunk(simulation; duration=5000),
    ]
end

# ! the first 0-v chunk is important to create the hex pattern

h = run_simulation(simulation, chunks; frame_every_n=20)
 