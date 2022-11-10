using GeneralAttractors
using Plots

using GeneralAttractors.Simulations
using GeneralAttractors.Kernels
using GeneralAttractors.ManifoldUtils
using GeneralAttractors: lerp
import GeneralAttractors.ManifoldUtils: sphere_embedding
using Distances

# --------------------------------- make net --------------------------------- #
n = (64, 64)
function ξ_s(i::Int, j::Int)::Vector
    [lerp(i, n[1], -π, π), lerp(j, n[2], -π / 2, π / 2)]
end
d_s = SphericalAngle()
# k_s = DiffOfExpKernel(; λ = 0.75)
k_s = MexicanHatKernel(; α = 0.1, σ = 0.5)

cover = CoverSpace(S², S², (x, y) -> [x, y])

can = CAN(
    "sphere",
    cover,
    n,
    ξ_s,
    d_s,
    k_s;
    offset_size = 1.15,
    φ = sphere_embedding,
    # offsets=O,
    # Ω=Ω
)



# --------------------------------- simulate --------------------------------- #
dt = 0.5
duration = 250  # ms   

nframes = (Int ∘ round)(duration / dt)
trajectory = Trajectory(can; T = nframes, σθ = 0.0, θ₀ = deg2rad(0), σv = 0.0, μv = 0.1)
simulation = Simulation(can, trajectory; b₀ = 1.0, η = 0.0)

# h = run_simulation(
#     simulation,
#     frame_every_n = 20,
#     discard_first_ms = 0,
#     average_over_ms = 1,
# )
