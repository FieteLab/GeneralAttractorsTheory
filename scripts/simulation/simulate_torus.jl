using GeneralAttractors
using Plots

using Distances
import GeneralAttractors: lerp
using GeneralAttractors.Kernels

@info "     ... torus attractor"
# neurons position and distance function
n = (64, 64)
function ξ_t(i::Int, j::Int)::Vector  # neurons coordinates function
    n̂_i, n̂_j = Int(n[1] / 2), Int(n[2] / 2)
    [lerp(i, n[1], -n̂_i, n̂_i), lerp(j, n[2], -n̂_j, n̂_j)]   # ∈ [-n/2, n/2] × [-n/2, n/2]
end
d_t = PeriodicEuclidean([n...])  # distance function over a torus manifold
# connectivity kernel

k_t = DiffOfExpKernel(; λ = 13.0)

# construct network
T = CAN("torus", n, ξ_t, d_t, k_t; λ=0.50)


@info "starting simulation"
simulation = Simulation(T; η = 0.1, b₀=1)


MODE = :CONSTANT  # random or constant

if MODE == :CONSTANT
    chunks =
        map(
            θ -> ConstantChunk([cos(θ), sin(θ)], simulation; duration = 250),
            [π + π/4, 2π],
            # [0, π / 4, π / 2, 3 / 4 * π, π, -3 / 4 * π, -π / 2, -π / 4],
        ) |> collect
    chunks = [RandomChunk(simulation; duration = 200, μ₀ = 1.0, σ = 1), chunks...]
else
    chunks = [RandomChunk(simulation; duration = 150_000, μ₀ = 1.0, σ = 1)]
end
# plot(chunks[1])

h = run_simulation(
    simulation,
    chunks;
    frame_every_n = MODE == :CONSTANT ? 20 : nothing,
    discard_first_ms = MODE == :CONSTANT ? 0 : 5000,
    average_over_ms = 20,
)
