using GeneralAttractors.Simulations

using GeneralAttractors.ManifoldUtils: Mobius, ψ_t, ψ_θ1
using Plots

x₀ = [0.3, 0]

t = Trajectory(
    Mobius();
    T = 500,
    σ = [1, 0, 1],
    x₀ = x₀,
    still = still,
    vmax = 0.03,
    modality = :constant,
    n_piecewise_segments = 3,
)

p = plot(
    # xlim=[-0.75, 0.75], ylim=[-0.5, 2π+0.5], 
    aspect_ratio = :equal,
    size = (400, 800),
    grid = false,
)

# plot traj
plot!(t.X[1:10:end, 1], t.X[1:10:end, 2], lw = 2, color = :black, label = "traj")


for t = -1/2:0.4:1/2, θ = 0:0.25:2π
    d1 = ψ_t(t, θ) .* 0.15
    plot!([t, t + d1[1]], [θ, θ + d1[2]], color = :black, label = nothing)
    plot!([t + 1.5, t + d1[1] + 1.5], [θ, θ + d1[2]], color = :black, label = nothing)

    d1 = ψ_θ1(t, θ) .* 0.15
    plot!([t, t + d1[1]], [θ, θ + d1[2]], color = :red, label = nothing)

    d1 = ψ_θ2(t, θ) .* 0.15
    plot!([t + 1.5, t + d1[1] + 1.5], [θ, θ + d1[2]], color = :green, label = nothing)
end

p
