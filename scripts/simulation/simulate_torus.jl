using Plots


using GeneralAttractors
using GeneralAttractors.Simulations


using Distances
using GeneralAttractors.Kernels
using GeneralAttractors: lerp
using GeneralAttractors.ManifoldUtils


# include("../networks/torus.jl")

# --------------------------------- simulate --------------------------------- #
dt = 0.5
duration = 300
still = 100  # initialization period        

# select neurons to initialize
x₀ = [0, 0]
d = map(i -> toruscan.metric(x₀, toruscan.X[:, i]), 1:size(toruscan.X, 2))
activate = zeros(length(d))
activate[d.<2] .= 1

# initialize trajectory and simulation
nframes = (Int ∘ round)(duration / dt)
trajectory = Trajectory(
        toruscan; T = nframes, dt=dt,
        σv = 0.01, μv = 0.05, vmax=0.1,
        σθ = 0.3, θ₀ = nothing, 
        still=still
)
simulation = Simulation(toruscan, trajectory; η = 0.5, b₀=0.5)

# run
h, X̄ = @time run_simulation(
    simulation;
    frame_every_n = 20,
    discard_first_ms = 0,
    average_over_ms = 1,
    fps = 10,
    s₀=1.0 .* activate,
);     

nothing

# TODO decoding has to be initialized at the end of the `still` phase with decoded peak position
# TODO fix mess in simulation plotting to be more consisten


