using GeneralAttractors
using Plots

using GeneralAttractors.Simulations
using GeneralAttractors.Kernels
using GeneralAttractors.ManifoldUtils
using GeneralAttractors: lerp
import GeneralAttractors.ManifoldUtils: sphere_embedding
using Distances
import LinearAlgebra: norm
import GeneralAttractors.Simulations: plot_trajectory_and_decoded



include("../networks/sphere.jl")

# --------------------------------- simulate --------------------------------- #
dt = 0.5
duration = 500  # ms   
x₀ = [1, 0, 0]
x₀ /= norm(x₀)
still = 100  # initialization period                                                                             
dmin = 0.5  # minimal distance from x₀ for state intialization


nframes = (Int ∘ round)(duration / dt)
trajectory = Trajectory(
    spherecan;
    T = nframes,
    x₀ = x₀,
    vmax = 0.01,
    still = still,
    σv = 0,
    μv = [0, 0, -1]
)

plot(trajectory) |> display
simulation = Simulation(spherecan, trajectory; b₀ = 0.25, η = 0.0)


# get activation to initialize bump
activate = map(p -> euclidean(x₀, p) < dmin, eachcol(spherecan.X)) .* 1

# simulate
h, X̄ = @time run_simulation(
    simulation,
    frame_every_n = 10,
    discard_first_ms = 0,
    average_over_ms = 10,
    fps = 10,
    s₀ = activate,
    savefolder = "sphere",
    savename = "decoding",
);


plot_trajectory_and_decoded(trajectory, X̄) |> display
