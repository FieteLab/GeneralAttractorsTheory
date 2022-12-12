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

include("../networks/mobius.jl")
println(mobiuscan.n)

# --------------------------------- simulate --------------------------------- #
dt = 0.5
duration = 200
still = 50  # initialization period  
dmin = 0.5  # minimal distance from x₀ for state intialization  

# select neurons to initialize
x₀ = [0.0, 5]
d = map(i -> mobiuscan.metric(x₀, mobiuscan.X[:, i]), 1:size(mobiuscan.X, 2))
activate = zeros(length(d))
activate[d.<dmin] .= 1

# initialize trajectory and simulation
nframes = (Int ∘ round)(duration / dt)
trajectory = Trajectory(
    mobiuscan;
    T = nframes,
    dt = dt,
    σv = [0, .5, .5],
    μv = [.1, 0, 0],
    x₀ = x₀,
    still = still,
    vmax = 0.1,
    scale = 1
)
simulation = Simulation(mobiuscan, trajectory; 
        η = 0.0, b₀ = 0.5, τ = 5.0)


# run
h, X̄ = @time run_simulation(
    simulation;
    frame_every_n = 10,
    discard_first_ms = 0,
    average_over_ms = 10,
    fps = 10,
    s₀ = 1.0 .* activate,
    φ = mobius_embedding,
    savefolder = "mobius",
    savename = "test",
);

# plot_trajectory_and_decoded(trajectory, X̄) |> display
nothing
