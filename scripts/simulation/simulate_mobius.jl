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
duration = 250
still = 50  # initialization period        

# select neurons to initialize
# x₀ = [0, 0]
# d = map(i -> mobiuscan.metric(x₀, mobiuscan.X[:, i]), 1:size(mobiuscan.X, 2))
# activate = zeros(length(d))
# activate[d.<2] .= 1

# initialize trajectory and simulation
nframes = (Int ∘ round)(duration / dt)
trajectory = Trajectory(
    mobiuscan; 
    T = nframes, 
    # x₀=x₀,
    vmax=0.0035,
    still=still,
    modality=:piecewise,
    n_piecewise_segments=2,
    σ=[1, 1, 1]
)
simulation = Simulation(mobiuscan, trajectory; η = 0.0, b₀=0.31)

# run
h, X̄ = @time run_simulation(
    simulation;
    frame_every_n = 20,
    discard_first_ms = 0,
    average_over_ms = 1,
    fps = 10,
    # s₀=1.0 .* activate,
);     

plot_trajectory_and_decoded(trajectory, X̄) |> display

nothing

# TODO fix visuals
# TODO fix decoding