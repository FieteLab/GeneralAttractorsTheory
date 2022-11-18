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
duration = 2000  # ms   
x₀ = [1, -1, 0]
x₀ /= norm(x₀)
still = 50  # initialization period                                                                             
dmin = 0.5  # minimal distance from x₀ for state intialization


nframes = (Int ∘ round)(duration / dt)
trajectory = Trajectory(
    spherecan; 
    T = nframes, 
    x₀=x₀,
    vmax=0.00375,
    still=still,
    modality=:piecewise,
    n_piecewise_segments=5,
    # σ=[1, 1, 1]
)

# get activation to initialize bump
activate = map(p -> euclidean(x₀, p) < dmin, eachcol(spherecan.X)) .* 1

# simulate
simulation = Simulation(spherecan, trajectory; b₀ = 0.27, η = 0.0, )
h, X̄ = @time run_simulation(
    simulation,
    frame_every_n = 20,
    discard_first_ms = 0,
    average_over_ms = 10,
    fps=10,
    s₀=activate,
    savefolder="sphere",
    savename="decoding",
);


plot_trajectory_and_decoded(trajectory, X̄) |> display
