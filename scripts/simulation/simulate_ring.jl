using Plots
using Term
install_term_stacktrace()

using GeneralAttractors
using GeneralAttractors.Simulations


using Distances
using GeneralAttractors.Kernels
using GeneralAttractors: lerp
using GeneralAttractors.ManifoldUtils
import GeneralAttractors.Simulations: plot_trajectory_and_decoded

include("../networks/ring.jl")

# --------------------------------- simulate --------------------------------- #
dt = 0.5
duration = 500
still = 50  # initialization period        


θ₀ = 3.14 # initialize state at center of mfld
d = ringcan.metric.(θ₀, ringcan.X[1, :])
activate = zeros(length(d))
activate[d.<0.1] .= 1

# initialize trajectory and simulation
nframes = (Int ∘ round)(duration / dt)
trajectory = Trajectory(
    ringcan;
    T = nframes,
    dt = dt,
    σθ = 0.0,
    θ̇₀ = 0.05,
    θ₀ = θ₀,
    still = still,
)
simulation = Simulation(ringcan, trajectory; η = 0.0, b₀ = 0.05)


h, X̄ = @time run_simulation(
    simulation;
    frame_every_n = 20,
    discard_first_ms = 0,
    average_over_ms = 10,
    fps = 10,
    s₀ = 1.0 .* activate,
    savefolder = "ring",
    savename = "test",
);

# # plot_trajectory_and_decoded(trajectory, X̄) |> display
nothing
