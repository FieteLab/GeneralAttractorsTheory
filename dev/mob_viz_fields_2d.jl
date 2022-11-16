using GeneralAttractors.Simulations

using GeneralAttractors.ManifoldUtils: Mobius, ψ_t, ψ_θ1
using Plots


t = Trajectory(Mobius(); T=500,     
    # modality=:constant, σ=[0, -.1, 0], 
    vmax=0.01)



p = plot(t.X[:, 1], t.X[:, 2], xlim=[-0.75, 0.75], ylim=[-0.5, 2π+0.5], lw=2, color=:black, label="traj", 
        aspect_ratio=:equal, size=(400, 800), grid=false)


for t in -1/2:.2:1/2, θ in 0:.5:2π
    d1 = ψ_t(t, θ) .* .15
    plot!([t, t+d1[1]], [θ, θ+d1[2]], color=:black, label=nothing)

    d1 = ψ_θ1(t, θ) .* .15
    plot!([t, t+d1[1]], [θ, θ+d1[2]], color=:red, label=nothing)

end

p