using Plots


using GeneralAttractors
using GeneralAttractors.Simulations


using Distances
using GeneralAttractors.Kernels
using GeneralAttractors: lerp
using GeneralAttractors.ManifoldUtils


include("../networks/torus.jl")

# --------------------------------- simulate --------------------------------- #
dt = 0.5
duration = 300

nframes = (Int ∘ round)(duration / dt)
trajectory = Trajectory(toruscan; T = nframes, μ = 0.25, σθ = 0.3, σv = 0.5, θ₀ = nothing)
simulation = Simulation(toruscan, trajectory; η = 0.5, b₀=0.5)

h, X̄ = @time run_simulation(
    simulation;
    frame_every_n = 20,
    discard_first_ms = 0,
    average_over_ms = 1,
    fps = 10,
);     

nothing

# TODO make decoding work for higher degree cover space.
# TODO fix mess in simulation plotting to be more consisten


