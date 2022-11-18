using Plots


using GeneralAttractors
using GeneralAttractors.Simulations

using Term
using Distances
using GeneralAttractors.Kernels
using GeneralAttractors: lerp
using GeneralAttractors.ManifoldUtils
import GeneralAttractors.Simulations: plot_trajectory_and_decoded

include("../scripts/networks/torus.jl")

# --------------------------------- simulate --------------------------------- #
dt = 0.5
duration = 3_000
still = 50  # initialization period        
n_sims = 200


x0, x1 = minimum(toruscan.X[2, :]), maximum(toruscan.X[2, :])
for i in 1:n_sims
    println(Panel("Starting condition $i/$n_sims", style="red", justify=:center))

    # select neurons to initialize
    x₀ = [rand(x0:1:x1), rand(x0:1:x1)]
    d = map(i -> toruscan.metric(x₀, toruscan.X[:, i]), 1:size(toruscan.X, 2))
    activate = zeros(length(d))
    activate[d.<2] .= 1

    # initialize trajectory and simulation
    nframes = (Int ∘ round)(duration / dt)
    trajectory = Trajectory(
            toruscan; T = nframes, dt=dt,
            σv = 0.00, μv = 0.05, vmax=0.05,
            σθ = 0.2, θ₀ = nothing, 
            still=still
    )
    simulation = Simulation(toruscan, trajectory; η = 0.0, b₀=0.31)

    # run
    h, X̄ = @time run_simulation(
        simulation;
        frame_every_n = nothing,
        discard_first_ms = 100,
        average_over_ms = 10,
        fps = 10,
        s₀=1.0 .* activate,
        savefolder="abstract",
        savename="torus_$i",
    );     

    # plot_trajectory_and_decoded(trajectory, X̄) |> display
end
