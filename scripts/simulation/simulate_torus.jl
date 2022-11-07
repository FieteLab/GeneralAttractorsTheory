using Plots


using GeneralAttractors
using GeneralAttractors.Simulations


using Distances
using GeneralAttractors.Kernels
using GeneralAttractors: lerp

n = (64, 64)
function ξ_t(i::Int, j::Int)::Vector  # neurons coordinates function
    n̂_i, n̂_j = Int(n[1] / 2), Int(n[2] / 2)
    [lerp(i, n[1], -n̂_i, n̂_i), lerp(j, n[2], -n̂_j, n̂_j)]   # ∈ [-n/2, n/2] × [-n/2, n/2]
end
d_t = PeriodicEuclidean([n...])  # distance function over a torus manifold
 
# connectivity kernel
k_t = DiffOfExpKernel(; λ = 13.0)
# k_t = LocalGlobalKernel(α = 0.1, σ = 20.0, β = 0.005)  

T = CAN("torus", ℝ², n, ξ_t, d_t, k_t;) 


trajectory = Trajectory(T; T = 1000, μ = 0.4, θ = 0.45)
simulation = Simulation(T, trajectory; η = 0.0)


h = @time run_simulation(
    simulation;
    frame_every_n = 20,
    discard_first_ms = 0,
    average_over_ms = 20,
    fps=10
)
