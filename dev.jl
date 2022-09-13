using GeneralAttractors
import GeneralAttractors.Simulations: AbstractChunk
simulation = Simulation(torus_attractor)
using Statistics

S = load_simulation_history("ring_sim").S


s1 = mean(reshape(S, (256, 2, 400, 5)); dims=4)