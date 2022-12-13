using Plots


using GeneralAttractors
using GeneralAttractors.Simulations
using Term
install_term_stacktrace()

using Distances
using GeneralAttractors.Kernels
using GeneralAttractors: lerp
using GeneralAttractors.ManifoldUtils
import GeneralAttractors.Simulations: plot_trajectory_and_decoded

include("../networks/torus.jl")

# --------------------------------- simulate --------------------------------- #
dt = 0.5
duration = 200
still = 100  # initialization period        

x₀ = [3.15, 3.15] # initialize state at center of mfld
d = map(i -> toruscan.metric(x₀, toruscan.X[:, i]), 1:size(toruscan.X, 2))
activate = zeros(length(d))
activate[d.<2.5] .= 1

# initialize trajectory and simulation
nframes = (Int ∘ round)(duration / dt)
trajectory = Trajectory(
    toruscan;
    T = nframes,
    dt = dt,
    σv = 0,
    μv = [0, 1],
    x₀ = x₀,
    vmax = 0.05,
    still = still,
    scale = 2,
)
simulation = Simulation(toruscan, trajectory; η = 0.0, b₀ = 0.5, τ = 5.0)

plot(trajectory)

# run
h, X̄ = @time run_simulation(
    simulation;
    frame_every_n = 20,
    discard_first_ms = 0,
    average_over_ms = 10,
    fps = 4,
    s₀ = 1.0 .* activate,
    savefolder = "torus",
    savename = "test",
);

# plot_trajectory_and_decoded(trajectory, X̄) |> display
nothing
