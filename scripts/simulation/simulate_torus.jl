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
duration = 800

nframes = (Int ∘ round)(duration / dt)
trajectory = Trajectory(toruscan; T = nframes, μ = 0.5, σθ = 0.3, σv = 0.5, θ₀ = nothing)
simulation = Simulation(toruscan, trajectory; η = 0.0)


h = @time run_simulation(
    simulation;
    frame_every_n = 20,
    discard_first_ms = 250,
    average_over_ms = 1,
    fps = 10,
)
