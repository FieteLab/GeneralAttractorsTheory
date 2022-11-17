using Plots


using GeneralAttractors
using GeneralAttractors.Simulations


using Distances
using GeneralAttractors.Kernels
using GeneralAttractors: lerp
using GeneralAttractors.ManifoldUtils
import GeneralAttractors.Simulations: plot_trajectory_and_decoded

# include("../networks/torus.jl")

# --------------------------------- simulate --------------------------------- #
dt = 0.5
duration = 2000
still = 50  # initialization period        


y0, y1 = minimum(toruscan.X[2, :]), maximum(toruscan.X[2, :])
for (i, y) in enumerate(range(y0, y1, length=20))
    # select neurons to initialize
    if i%2 == 0
        x₀ = [0, y]
        θ   = 0
    else
        x₀ = [y, 0]
        θ   = π/2
    end
    d = map(i -> toruscan.metric(x₀, toruscan.X[:, i]), 1:size(toruscan.X, 2))
    activate = zeros(length(d))
    activate[d.<2] .= 1

    # initialize trajectory and simulation
    nframes = (Int ∘ round)(duration / dt)
    trajectory = Trajectory(
            toruscan; T = nframes, dt=dt,
            σv = 0.00, μv = 0.05, vmax=0.05,
            σθ = 0, θ₀ = θ, 
            still=still
    )
    simulation = Simulation(toruscan, trajectory; η = 0.0, b₀=0.31)

    # run
    h, X̄ = @time run_simulation(
        simulation;
        frame_every_n = 10,
        discard_first_ms = 100,
        average_over_ms = 10,
        fps = 10,
        s₀=1.0 .* activate,
        savefolder="abstract",
        savename="torus_$i",
    );     

    # plot_trajectory_and_decoded(trajectory, X̄) |> display
    nothing
end
