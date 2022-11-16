using GeneralAttractors.Simulations

using GeneralAttractors.ManifoldUtils: Mobius, ψ_t, ψ_θ1
using Plots


t = Trajectory(Mobius(); T=500,  
    σθ=0.0,   
    θ₀ = π/4,
    μv = 0.02,
    vmax=0.02)


p = plot(
    # xlim=[-0.75, 0.75], ylim=[-0.5, 2π+0.5], 
    aspect_ratio=:equal, size=(400, 800), grid=false
)

# plot traj
plot!(t.X[1:10:end, 1], t.X[1:10:end, 2], lw=2, color=:black, label="traj")


for t in -1/2:.2:1/2, θ in 0:.25:2π
    d1 = ψ_t(t, θ) .* .05
    plot!([t, t+d1[1]], [θ, θ+d1[2]], color=:black, label=nothing)

    # d1 = ψ_θ1(t, θ) .* .15
    # plot!([t, t+d1[1]], [θ, θ+d1[2]], color=:red, label=nothing)
    
    d1 = -ψ_θ2(t, θ) .* .15
    plot!([t, t+d1[1]], [θ, θ+d1[2]], color=:green, label=nothing)
end

p