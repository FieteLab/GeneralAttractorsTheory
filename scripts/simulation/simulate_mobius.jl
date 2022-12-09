using Plots
using GeneralAttractors
using GeneralAttractors.Simulations
using Term
install_term_stacktrace()
install_term_logger()

using Distances
using GeneralAttractors.Kernels
using GeneralAttractors: lerp
using GeneralAttractors.ManifoldUtils
import GeneralAttractors.Simulations: plot_trajectory_and_decoded

# include("../networks/mobius.jl")

# --------------------------------- simulate --------------------------------- #
dt = 0.5
duration = 200
still = 100  # initialization period  
dmin = 0.25  # minimal distance from x₀ for state intialization  

# select neurons to initialize
x₀ = [0.3, 0]
d = map(i -> mobiuscan.metric(x₀, mobiuscan.X[:, i]), 1:size(mobiuscan.X, 2))
activate = zeros(length(d))
activate[d.<dmin] .= 1

# initialize trajectory and simulation
nframes = (Int ∘ round)(duration / dt)
trajectory = Trajectory(
    mobiuscan;
    T = nframes,
    σ = [1, 0.1, 0.1],
    x₀ = x₀,
    still = still,
    vmax = 0.0075,
    modality = :piecewise,
    n_piecewise_segments = 12,
)
simulation = Simulation(mobiuscan, trajectory; η = 0.0, b₀ = 0.20)

# run
h, X̄ = @time run_simulation(
    simulation;
    frame_every_n = 20,
    discard_first_ms = 0,
    average_over_ms = 10,
    fps = 10,
    s₀ = 1.0 .* activate,
    φ = mobius_embedding,
    savefolder = "mobius",
    savename = "test",
);

plot_trajectory_and_decoded(trajectory, X̄) |> display
nothing
