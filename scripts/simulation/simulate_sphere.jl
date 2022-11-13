using GeneralAttractors
using Plots

using GeneralAttractors.Simulations
using GeneralAttractors.Kernels
using GeneralAttractors.ManifoldUtils
using GeneralAttractors: lerp
import GeneralAttractors.ManifoldUtils: sphere_embedding
using Distances
import LinearAlgebra: norm

include("../networks/sphere.jl") 


# --------------------------------- simulate --------------------------------- #
dt = 0.5
duration = 600  # ms   
x₀ = [0, 1, 0]
x₀ /= norm(x₀)
still = 50  # initialization period
dmin = 0.5


nframes = (Int ∘ round)(duration / dt)
trajectory = Trajectory(spherecan; T = nframes, σ=[0.0, 0.0, 0.0], x₀=x₀, scale=0.05, still=still)


# get activation to initialize bump
activate = map(p -> euclidean(x₀, p) < dmin, eachcol(spherecan.X)) .* 1
sum(activate/length(activate)) |> println

# simulate
simulation = Simulation(spherecan, trajectory; b₀ = 0.1, η = 0.0, )
h, X̄ = run_simulation(
    simulation,
    frame_every_n = 20,
    discard_first_ms = 0,
    average_over_ms = 2,
    fps=5,
    s₀=1.0 .* activate,
    activation_steps=still
);
