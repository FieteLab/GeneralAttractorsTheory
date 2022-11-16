using Plots


using GeneralAttractors
using GeneralAttractors.Simulations


using Distances
using GeneralAttractors.Kernels
using GeneralAttractors: lerp
using GeneralAttractors.ManifoldUtils
import GeneralAttractors.Simulations: plot_trajectory_and_decoded

include("../networks/mobius.jl")

# --------------------------------- simulate --------------------------------- #
dt = 0.5
duration = 400
still = 100  # initialization period  
dmin = 0.25  # minimal distance from x₀ for state intialization  

# select neurons to initialize
x₀ = [-0.3, .5]
d = map(i -> mobiuscan.metric(x₀, mobiuscan.X[:, i]), 1:size(mobiuscan.X, 2))
activate = zeros(length(d))
activate[d.<dmin] .= 1

# initialize trajectory and simulation
nframes = (Int ∘ round)(duration / dt)
trajectory = Trajectory(
    mobiuscan; 
    T = nframes, 
    σ = [1, 0, 0],
    x₀=x₀,
    still=still,
    vmax=0.0075,
    modality=:constant,
    n_piecewise_segments=3,
)
simulation = Simulation(mobiuscan, trajectory; η = 0.0, b₀=0.3)

# run
h, X̄ = @time run_simulation(
    simulation;
    frame_every_n = 20,
    discard_first_ms = 0,
    average_over_ms = 1,
    fps = 10,
    s₀=1.0 .* activate,
    φ=mobius_embedding,
);     

# plot_trajectory_and_decoded(trajectory, X̄) |> display

nothing

# TODO fix visuals
# TODO fix decoding