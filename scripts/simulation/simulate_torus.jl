 using Plots


using GeneralAttractors
using GeneralAttractors.Simulations
using Term
install_term_stacktrace(hide_frames = true)

using Distances
using GeneralAttractors.Kernels
using GeneralAttractors: lerp, moving_average
using GeneralAttractors.ManifoldUtils
import GeneralAttractors.Simulations: plot_trajectory_and_decoded, plot_on_mfld_trajectory_and_history
import GeneralAttractors.Analysis: get_bump_speed

include("../networks/torus.jl")

# --------------------------------- simulate --------------------------------- #
dt = 0.5
duration = 2500
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
    σv = .15,
    # σv = 0,
    # μv = [.1, .001],
    x₀ = x₀,
    vmax = 0.05,
    still = still,
    scale = 12,
    smoothing_window=1000,
    # modality=:piecewise,
    # n_piecewise_segments=10
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
traj_speed = norm.(eachcol(h.v))[20:end] / h.Δt
mfld_speed = moving_average(get_bump_speed(toruscan, h)[1], 31)

p3 = plot(cumsum(traj_speed), lw=4, color=:black, label="trajectory", title="distance travelled")
plot!(p3, cumsum(mfld_speed), lw=2, color=:red, label="mfld", alpha=.8)



p4 = plot(traj_speed, lw=4, color=:black, label="input trajectory", title="speed")
# plot!(p4, mfld_speed, lw=0.5, color=:red, label="mfld speed", alpha=.25)
plot!(p4, moving_average(mfld_speed, 25), lw=2, color=:red, label="mfld", alpha=.8)


p5 = scatter(traj_speed, mfld_speed, 
    label=nothing, color=:black, 
    aspect_ratio=:equal, xlabel="trajectory", ylabel="mfld", title="input vs bump speed"); 
plot!([0, .05], [0, .05], lw=3, alpha=.5, ls=:dash, color=:black, label=nothing)

plot(p1, p3, p4,p5,  layout=(4, 2), size=(800, 800)) |> display
nothing
