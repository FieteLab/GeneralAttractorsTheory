using Plots


using GeneralAttractors
using GeneralAttractors.Simulations


using Distances
using GeneralAttractors.Kernels
using GeneralAttractors: lerp
using GeneralAttractors.ManifoldUtils
import GeneralAttractors.Simulations: plot_trajectory_and_decoded

include("../scripts/networks/mobius.jl")

# --------------------------------- simulate --------------------------------- #
dt = 0.5
duration = 3000
still = 50  # initialization period  
dmin = 0.25  # minimal distance from x₀ for state intialization  

n_sims = 100

for i in 1:n_sims
    println(Panel("Starting condition $i/$n_sims", style="red", justify=:center))

    # select neurons to initialize
    x₀ = [rand(-1/2:.1:1/2), rand(0:.1:2π)]
    d = map(i -> mobiuscan.metric(x₀, mobiuscan.X[:, i]), 1:size(mobiuscan.X, 2))
    activate = zeros(length(d))
    activate[d.<dmin] .= 1

    # initialize trajectory and simulation
    nframes = (Int ∘ round)(duration / dt)
    trajectory = Trajectory(
        mobiuscan; 
        T = nframes, 
        σ = [1, 1, 1],
        x₀=x₀,
        still=still,
        vmax=0.0077,
        modality=:piecewise,
        n_piecewise_segments=3, 
    )
    simulation = Simulation(mobiuscan, trajectory; η = 0.0, b₀=0.3)

    # run
    h, X̄ = @time run_simulation(
        simulation;
        frame_every_n = nothing,
        discard_first_ms = 0,
        average_over_ms = 10,
        fps = 10,
        s₀=1.0 .* activate,
        φ=mobius_embedding,
        savefolder="abstract",
        savename="mobius_$(i+60)",
    );     
end