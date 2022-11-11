using GeneralAttractors
using Plots

using GeneralAttractors.Simulations
using GeneralAttractors.Kernels
using GeneralAttractors.ManifoldUtils
using GeneralAttractors: lerp
import GeneralAttractors.ManifoldUtils: sphere_embedding
using Distances


include("../networks/sphere.jl") 


# --------------------------------- simulate --------------------------------- #
dt = 0.5
duration = 1500  # ms   

nframes = (Int ∘ round)(duration / dt)
trajectory = Trajectory(spherecan; T = nframes, σ=[0.0, 0.0, 0.0], x₀=[-π/2, 0], still=300)
simulation = Simulation(spherecan, trajectory; b₀ = 0.1, η = 0.0)


plot(trajectory.X[1:end, 1])

h = run_simulation(
    simulation,
    frame_every_n = 20,
    discard_first_ms = 0,
    average_over_ms = 2,
    fps=5,
)
