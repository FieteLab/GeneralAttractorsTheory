using Plots
using GeneralAttractors
using GeneralAttractors.Simulations
using Term

using Distances
using GeneralAttractors.Kernels
using GeneralAttractors: lerp
using GeneralAttractors.ManifoldUtils
import GeneralAttractors.Simulations: plot_trajectory_and_decoded


can = mobius_maker(:default)

# --------------------------------- simulate --------------------------------- #
dt = 0.5
duration = 400
still = 50  # initialization period  
dmin = 0.5  # minimal distance from x₀ for state intialization  

# select neurons to initialize
x₀ = [0.5, 5]
d = map(i -> can.metric(x₀, can.X[:, i]), 1:size(can.X, 2))
activate = zeros(length(d))
activate[d.<dmin] .= 1

# initialize trajectory and simulation
nframes = (Int ∘ round)(duration / dt)
trajectory = Trajectory(
    can;
    T = nframes,
    dt = dt,
    σv = [0.5, 0.3],
    # σv = 0,
    μv = [0.1, 0.0],
    x₀ = x₀,
    still = still,
    vmax = 0.1,
)
plot(trajectory) |> display


simulation = Simulation(can, trajectory; η = 0.0, b₀ = 0.5, τ = 5.0)



# run
h, X̄ = @time run_simulation(
    simulation;
    # frame_every_n = 25,
    discard_first_ms = 0,
    average_over_ms = 10,
    s₀ = 1.0 .* activate,
    # φ = mobius_embedding,
    # savefolder = "mobius",
    # savename = "test",
);

plot_trajectory_and_decoded(trajectory, X̄) |> display
nothing
