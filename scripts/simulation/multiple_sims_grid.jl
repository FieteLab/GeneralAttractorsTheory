using GeneralAttractors
using Plots
using Term

using GeneralAttractors.Simulations
using GeneralAttractors.Kernels
using GeneralAttractors.ManifoldUtils
using GeneralAttractors: lerp
import GeneralAttractors.ManifoldUtils: sphere_embedding
using Distances
import LinearAlgebra: norm
import GeneralAttractors.Simulations: plot_trajectory_and_decoded
import GeneralAttractors.ManifoldUtils: fibonacci_sphere


include("../scripts/networks/sphere.jl")

# --------------------------------- simulate --------------------------------- #
dt = 0.5
duration = 3_000  # ms  
still = 50  # initialization period                                                                             
dmin = 0.5  # minimal distance from x₀ for state intialization
nframes = (Int ∘ round)(duration / dt)

n_sims = 100


pts = fibonacci_sphere(n_sims)
for i = 1:n_sims
    println(Panel("Starting condition $i/$n_sims", style = "red", justify = :center))

    x₀ = pts[:, i]
    x₀ /= norm(x₀)

    trajectory = Trajectory(
        spherecan;
        T = nframes,
        x₀ = x₀,
        vmax = 0.0035,
        still = still,
        modality = :piecewise,
        n_piecewise_segments = 6,
    )

    # get activation to initialize bump
    activate = map(p -> euclidean(x₀, p) < dmin, eachcol(spherecan.X)) .* 1

    # simulate
    simulation = Simulation(spherecan, trajectory; b₀ = 0.20, η = 0.0)
    h, X̄ = @time run_simulation(
        simulation,
        frame_every_n = nothing,
        discard_first_ms = 0,
        average_over_ms = 10,
        fps = 10,
        s₀ = activate,
        savefolder = "abstract",
        savename = "sphere_$(i+50)",
    )

end
