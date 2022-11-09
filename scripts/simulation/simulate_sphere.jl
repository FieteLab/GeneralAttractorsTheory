using GeneralAttractors
using Plots

using GeneralAttractors.Simulations
using GeneralAttractors.Kernels
using GeneralAttractors.Manifolds
using GeneralAttractors: lerp
import GeneralAttractors.Manifolds: sphere_embedding
using Distances

# --------------------------------- make net --------------------------------- #
n = (64, 64)
function ξ_s(i::Int, j::Int)::Vector
    [lerp(i, n[1], -π, π), lerp(j, n[2], -π / 2, π / 2)]
end
d_s = SphericalAngle()
k_s = DiffOfExpKernel(; λ = 0.75)

cover = CoverSpace(S², S², (x, y) -> [mod(x, 1), mod(y, 1)])
sph = CAN("sphere", cover, n, ξ_s, d_s, k_s; offset_size=[0.2, 0.2, 0.1, 0.1], φ=sphere_embedding)



# --------------------------------- simulate --------------------------------- #
dt = 0.5
duration = 800

nframes = (Int ∘ round)(duration / dt)
trajectory = Trajectory(sph; T = nframes, )
simulation = Simulation(sph, trajectory)

h = run_simulation(
    simulation,
    frame_every_n = 20,
    discard_first_ms = 0,
    average_over_ms = 20,
)
