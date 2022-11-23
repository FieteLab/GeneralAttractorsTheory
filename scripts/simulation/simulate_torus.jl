using Plots


using GeneralAttractors
using GeneralAttractors.Simulations


using Distances
using GeneralAttractors.Kernels
using GeneralAttractors: lerp
using GeneralAttractors.ManifoldUtils
import GeneralAttractors.Simulations: plot_trajectory_and_decoded

include("../networks/torus.jl")

# --------------------------------- simulate --------------------------------- #
dt = 0.5
duration = 1000
still = 50  # initialization period        


x₀ = [0, 0] # inittialize state at center of mfld
d = map(i -> toruscan.metric(x₀, toruscan.X[:, i]), 1:size(toruscan.X, 2))
activate = zeros(length(d))
activate[d.<2] .= 1

# initialize trajectory and simulation
nframes = (Int ∘ round)(duration / dt)
trajectory = Trajectory(
    toruscan;
    T = nframes,
    dt = dt,
    σv = 0.00,
    μv = 0.1,
    vmax = 0.1,
    σθ = 0.2,
    θ₀ = nothing,
    still = still,
)
simulation = Simulation(toruscan, trajectory; η = 0.0, b₀ = 0.31)

# run
h, X̄ = @time run_simulation(
    simulation;
    frame_every_n = 10,
    discard_first_ms = 100,
    average_over_ms = 10,
    fps = 10,
    s₀ = 1.0 .* activate,
    savefolder = "torus",
    savename = "decoding",
);

plot_trajectory_and_decoded(trajectory, X̄) |> display
