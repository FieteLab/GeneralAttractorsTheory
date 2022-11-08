using Plots


using GeneralAttractors
using GeneralAttractors.Simulations


using Distances
using GeneralAttractors.Kernels
using GeneralAttractors: lerp
using GeneralAttractors.Manifolds


# ---------------------------------- network --------------------------------- #
cover = CoverSpace(ℝ², T, (x, y) -> [mod(x, 64), mod(y, 64)])


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
    OneForm(1, x -> 2sin(x/n[1])+0.25),
    OneForm(1, x -> -(2sin(x/n[1])+0.25)),
    OneForm(2, x -> 2sin(x/n[2])+0.25),
    OneForm(2, x -> -(2sin(x/n[2])+0.25)),
]

# make network
tor = CAN("torus", cover, n, ξ_t, d_t, k_t; Ω=Ω) 

# --------------------------------- simulate --------------------------------- #
dt= 0.5
duration = 800

nframes = (Int ∘ round)(duration/dt)
trajectory = Trajectory(tor; T = nframes, μ = 0.5, σθ = 0.0, σv = 0.0, θ₀=π/6)
simulation = Simulation(tor, trajectory; η = 0.0)


h = @time run_simulation(
    simulation;
    frame_every_n = 20,
    discard_first_ms = 250,
    average_over_ms = 20,
    fps=10
)
