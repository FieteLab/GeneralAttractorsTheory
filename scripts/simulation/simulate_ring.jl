using Plots


using GeneralAttractors
using GeneralAttractors.Simulations
using Term
Term.STACKTRACE_HIDDEN_MODULES[] = ["Plots"]
install_term_stacktrace()

using Distances
using GeneralAttractors.Kernels
using GeneralAttractors: lerp
using GeneralAttractors.ManifoldUtils
import GeneralAttractors.Simulations: plot_trajectory_and_decoded

include("../networks/ring.jl")

# --------------------------------- simulate --------------------------------- #
dt = 0.5
duration = 3000
still = 50  # initialization period        


θ₀ = 5 # initialize state at position
d = ringcan.metric.(θ₀, ringcan.X[1, :])
activate = zeros(length(d))
activate[d.<0.4] .= 1

# initialize trajectory and simulation
nframes = (Int ∘ round)(duration / dt)
trajectory = Trajectory(
    ringcan;
    T = nframes,
    dt = dt,
    σv = 0.5,
    μv = 0.0,
    x₀ = θ₀,
    still = still,
    vmax = 0.3,
    # scale=2,
)
simulation = Simulation(ringcan, trajectory; η = 0, b₀ = 0.1, τ = 10.0)


h, X̄ = @time run_simulation(
    simulation;
    frame_every_n = nothing,
    discard_first_ms = 0,
    average_over_ms = 10,
    fps = 10,
    s₀ = 1.0 .* activate,
    savefolder = "ring",
    savename = "test",
);

plot_trajectory_and_decoded(trajectory, X̄) |> display
nothing
