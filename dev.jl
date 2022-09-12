using GeneralAttractors
import GeneralAttractors.Can: CAN
using GeneralAttractors.Kernels
using Distances

# TODO test with non-square A
# TODO test with 1D networks


# neurons position and distance function
n = (64, 64)  
function ξ_t(i::Int, j::Int)::Vector  # neurons coordinates function
    p_i, p_j = (i-1)/(n[1]-1), (j-1)/(n[2]-1) # ∈ [0, 1]
    n̂_i, n̂_j = Int(n[1]/2), Int(n[2]/2)
    [
        -n̂_i*(1-p_i)+n̂_i*p_i, 
        -n̂_j*(1-p_j)+n̂_j*p_j, 
    ]  # ∈ [-n/2, n/2] × [-n/2, n/2]
end
d_t = PeriodicEuclidean([n...])  # distance function over a torus manifold

# connectivity kernel
k_t = DiffOfExpKernel(; λ = 13.0)

# construct network
can = CAN(n, ξ_t, d_t, k_t; offset_strength=1)


# chunks = [
#     SimulationChunk([cos(0), sin(0)]),
#     SimulationChunk([cos(π/5), sin(π/5)]),
#     SimulationChunk([cos(π/2-π/5), sin(π/2-π/5)]),
# ]

chunks = map(
    θ -> SimulationChunk([cos(θ), sin(θ)]; duration=600),
    [0, π/4, π/2, 3/4*π, π, -3/4*π, -π/2, -π/4]
) |> collect
chunks = SimulationChunk[SimulationChunk([0.0, 0.0]; duration=600), chunks...]

run_simulation(can, chunks)
