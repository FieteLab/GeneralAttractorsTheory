using GeneralAttractors
import GeneralAttractors.Simulations: AbstractChunk
simulation = Simulation(torus_attractor)


chunks = [
    RandomChunk(simulation)
]
run_simulation(simulation, chunks)
# plot(hcat(chunks[1].v...)[1, :])