using GeneralAttractors
using Plots

using Distances
import GeneralAttractors: lerp
using GeneralAttractors.Kernels


n = (64, 64)
function ξ_t(i::Int, j::Int)::Vector  # neurons coordinates function
    n̂_i, n̂_j = Int(n[1] / 2), Int(n[2] / 2)
    [lerp(i, n[1], -n̂_i, n̂_i), lerp(j, n[2], -n̂_j, n̂_j)]   # ∈ [-n/2, n/2] × [-n/2, n/2]
end
d_t = PeriodicEuclidean([n...])  # distance function over a torus manifold
k_t = DiffOfExpKernel(; λ = 13.0)

# construct network
T = CAN("torus", n, ξ_t, d_t, k_t; λ=1.0, σ=:relu)


@info "starting simulation"
simulation = Simulation(T; η = 0.1, b₀=1)


MODE = :RANDOM  # random or constant
if MODE == :CONSTANT
    # chunks =
    #     map(
    #         θ -> ConstantChunk([cos(θ), sin(θ)], simulation; duration = 150),
    #         [0, ], #π + π/4, 0, 2π],
    #         # [0, π / 4, π / 2, 3 / 4 * π, π, -3 / 4 * π, -π / 2, -π / 4],
    #     ) |> collect
    # chunks = [RandomChunk(simulation; duration = 200, μ₀ = 1.0, σ = 1), chunks...]

    chunks = [
        ConstantChunk([0, 0], simulation; duration=500),
        ConstantChunk([1.0, 0.0], simulation; duration=500),
        # ConstantChunk([0, 0], simulation; durat ion=250),
        # ConstantChunk([0, 0.8], simulation; duration=500),
        # ConstantChunk([0, 0], simulation; duration=250),
        # ConstantChunk([0.7, 0.7], simulation; duration=500),
    ]
else
    chunks = [
        ConstantChunk([0, 0], simulation; duration=250),
        ConstantChunk([0.5, 0], simulation; duration=400),
        ConstantChunk([0, 0.5], simulation; duration=400),
        ConstantChunk([0.5sin(π/4), 0.5cos(π/4)], simulation; duration=400),
        # RandomChunk(simulation; duration = 500, μ₀ = 10.5, σ = 1)
        ]
end

h = @time run_simulation(
    simulation,
    chunks;
    frame_every_n = MODE == :CONSTANT ? 20 : 20,
    discard_first_ms = MODE == :CONSTANT ? 0 : 0,
    average_over_ms = 20,
)
