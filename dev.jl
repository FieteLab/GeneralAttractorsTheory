using GeneralAttractors
import GeneralAttractors.Can: CAN
using GeneralAttractors.Kernels
using Distances


n = (48, 48)
function ξ_m(i::Int, j::Int)
    p_i, p_j = (i-1)/(n[1]-1), (j-1)/(n[2]-1) # ∈ [0, 1]
    [
        2π*p_i, 
        2π*p_j
    ]  # ∈ [0, 2π] × [0, 2π]
end
d_m = MobiusEuclidean(2π)

# connectivity kernel
k_m = DiffOfExpKernel(; λ = 13.0)

can = CAN(n, ξ_m, d_m, k_m; offset_size=0.1)
show_connectivity(can)




simulation = Simulation(
    can; frame_every_n=50, 
)

chunks = SimulationChunk[
    SimulationChunk([0.0, 0.0], simulation; duration=250), 
    SimulationChunk([1.0, 0.0], simulation; duration=200), 
    SimulationChunk([0.0, 0.0], simulation; duration=100), 
    SimulationChunk([0.0, 1.0], simulation; duration=200), 
]

run_simulation(simulation, chunks)
