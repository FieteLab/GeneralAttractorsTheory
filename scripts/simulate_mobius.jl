using GeneralAttractors
import GeneralAttractors.Can: CAN
using GeneralAttractors.Kernels
using Distances


simulation = Simulation(
    mobius_attractor;  
)

chunks = ConstantChunk[
    ConstantChunk([0.0, 0.0], simulation; duration=250), 
    ConstantChunk([1.0, 0.0], simulation; duration=800), 
    ConstantChunk([0.0, 0.0], simulation; duration=50), 
    ConstantChunk([0.0, 1.0], simulation; duration=800), 
    ConstantChunk([0.0, 0.0], simulation; duration=50), 
    ConstantChunk([-1.0, 0.0], simulation; duration=800), 
    ConstantChunk([0.0, 0.0], simulation; duration=50), 
    ConstantChunk([0.0, -1.0], simulation; duration=800), 
]

run_simulation(simulation, chunks; frame_every_n=10)
