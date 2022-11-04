using GeneralAttractors

using GeneralAttractors.Kernels
import GeneralAttractors: lerp
using Distances

@info "     ... ring attractor"
# neurons position and distance function
n = (256,)  # number of neurons in the ring
ξ_r(i::Int)::Vector = [lerp(i, n[1], 0.0, 2π)]  # neurons coordinates function
d_r = PeriodicEuclidean([2π])  # distance function

# construct network
attr = CAN("ring", n, ξ_r, d_r, MexicanHatKernel(α = 0.01, σ = .1); λ=1.5)


simulation = Simulation(attr; η = 0.0, b₀ = 1.0)


# chunks =
#     map(
#         θ -> ConstantChunk([Float64(θ)], simulation; duration = 250), 
#         [1]
#     ) |> collect

chunks = [
    ConstantChunk([0.0], simulation; duration = 500),
    ConstantChunk([3.5], simulation; duration = 750),
    ConstantChunk([-0.5], simulation; duration = 750),
]

run_simulation(simulation, chunks; frame_every_n = 20)
