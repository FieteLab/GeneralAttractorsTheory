using Plots


using GeneralAttractors
using GeneralAttractors.Simulations


using Distances
using GeneralAttractors.Kernels
using GeneralAttractors: lerp
using GeneralAttractors.Manifolds


# ---------------------------------- network --------------------------------- #
cover = CoverSpace(ℝ², T, (x, y) -> [mod(x - 32, 64), mod(y - 32, 64)])


n = (64, 64)
function ξ_t(i::Int, j::Int)::Vector  # neurons coordinates function
    n̂_i, n̂_j = Int(n[1] / 2), Int(n[2] / 2)
    [lerp(i, n[1], -n̂_i, n̂_i), lerp(j, n[2], -n̂_j, n̂_j)]   # ∈ [-n/2, n/2] × [-n/2, n/2]
end
d_t = PeriodicEuclidean([n...])  # distance function over a torus manifold

# connectivity kernel 
k_t = DiffOfExpKernel(; λ = 13.0)


# one forms
Ω = OneForm[
    OneForm(1, x -> sin(x / n[1]) + 1.25),
    OneForm(1, x -> -(sin(x / n[1]) + 1.25)),
    OneForm(2, x -> sin(x / n[2]) + 1.25),
    OneForm(2, x -> -(sin(x / n[2]) + 1.25)),
]

# make network
tor = CAN("torus", cover, n, ξ_t, d_t, k_t; Ω = Ω)

# --------------------------------- simulate --------------------------------- #
dt = 0.5
duration = 400

nframes = (Int ∘ round)(duration / dt)
trajectory = Trajectory(tor; T = nframes, μ = 0.5, σθ = 0.5, σv = 0.5, θ₀ = nothing)
simulation = Simulation(tor, trajectory; η = 0.0)


h = @time run_simulation(
    simulation;
    frame_every_n = 20,
    discard_first_ms = 50,
    average_over_ms = 20,
    fps = 10,
)


# plot(simulation, duration, nframes, 
#         trajectory.X[end, :], trajectory.V[end, :]; 
#         show_one_forms=true, dx=10, scale=2)
