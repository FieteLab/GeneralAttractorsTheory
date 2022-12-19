 using Plots


using GeneralAttractors
using GeneralAttractors.Simulations
using Term
install_term_stacktrace(hide_frames=true)

using Distances
using GeneralAttractors.Kernels
using GeneralAttractors: lerp, moving_average
using GeneralAttractors.ManifoldUtils
import GeneralAttractors.Simulations: plot_trajectory_and_decoded, plot_on_mfld_trajectory_and_history
import GeneralAttractors.Analysis: get_bump_speed

include("../networks/torus.jl")

# --------------------------------- simulate --------------------------------- #
dt = 0.5
duration = 1000
still = 100  # initialization period        

x₀ = [3.15, 3.15] # initialize state at center of mfld
d = map(i -> toruscan.metric(x₀, toruscan.X[:, i]), 1:size(toruscan.X, 2))
activate = zeros(length(d))
activate[d.<1.5] .= 1

# initialize trajectory and simulation
nframes = (Int ∘ round)(duration / dt)
trajectory = Trajectory(
    toruscan;
    T = nframes,
    dt = dt,
    σv = .1,
    # μv = [.1, .0],
    x₀ = x₀,
    vmax = 0.05,
    still = still,
    scale = 12,
    smoothing_window=100,
)
plot(trajectory) |> display

simulation = Simulation(toruscan, trajectory; η = 0.0, b₀ = 0.5, τ = 5.0)


# run
h, X̄ = @time run_simulation(
    simulation;
    frame_every_n = nothing,
    discard_first_ms = 0,
    average_over_ms = 1,
    fps = 4,
    s₀ = 1.0 .* activate,
    savefolder = "torus",
    savename = "test",
);

p1 = plot_trajectory_and_decoded(trajectory, X̄)
p2 = plot_on_mfld_trajectory_and_history(toruscan, trajectory, h)

# plot input vs bump speed
traj_speed = norm.(eachcol(h.v))[20:end]
mfld_speed, _ = get_bump_speed(toruscan, h)

p3 = plot(traj_speed, lw=4, color=:black, label="input trajectory")
plot!(p3, mfld_speed, lw=2, color=:red, label="mfld speed", alpha=.5)
plot!(p3, moving_average(mfld_speed, 11), lw=2, color=:green, label="mfld speed", alpha=.5)


plot(p1, p2, p3, layout=(3, 1), size=(800, 800)) |> display
nothing
