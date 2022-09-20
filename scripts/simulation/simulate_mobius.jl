using GeneralAttractors
import GeneralAttractors.Can: CAN
using GeneralAttractors.Kernels

# using Distances
simulation = Simulation(
    mobius_attractor;   b₀=1.0
)

chunks = ConstantChunk[
    ConstantChunk([0.0, 0.0], simulation; duration=250), 
    ConstantChunk([1.0, 0.0], simulation; duration=1000), 
    ConstantChunk([0.0, 1.0], simulation; duration=1000), 
    ConstantChunk([-1.0, 0.0], simulation; duration=1000), 
    ConstantChunk([0.0, -1.0], simulation; duration=1000), 
]

@time run_simulation(simulation, chunks; frame_every_n=25)


