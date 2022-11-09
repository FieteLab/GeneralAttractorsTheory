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

cover = CoverSpace(S², S², (x, y) -> [x, y])

O = [
    [cos(0), sin(0)],
    [cos(π * 1/4), sin(π * 1/4)],
    [cos(π * 2/4), sin(π * 2/4)],
    [cos(π * 3/4), sin(π * 3/4)],
    [cos(π), sin(π)],
    [cos(π+π*1/4), sin(π+π*1/4)],
    [cos(π+π*2/4), sin(π+π*2/4)],
    [cos(π+π*3/4), sin(π+π*3/4)],
]

Ω = OneForm[
    OneForm(1, (x, y) -> O[1]),
    OneForm(1, (x, y) -> O[2]),
    OneForm(1, (x, y) -> O[3]),
    OneForm(1, (x, y) -> O[4]),
    OneForm(2, (x, y) -> O[5]),
    OneForm(2, (x, y) -> O[6]),
    OneForm(2, (x, y) -> O[7]),
    OneForm(2, (x, y) -> O[8]),
]

# Ω = OneForm[
#     OneForm(1, (x, y) -> [1, 0]),
#     OneForm(1, (x, y) -> [1, 0]),
#     OneForm(1, (x, y) -> [0, 1]),
#     OneForm(1, (x, y) -> [0, 1]),
#     OneForm(2, (x, y) -> [-1, 0]),
#     OneForm(2, (x, y) -> [-1, 0]),
#     OneForm(2, (x, y) -> [0, -1]),
#     OneForm(2, (x, y) -> [0, -1]),
# ]


can = CAN("sphere", cover, n, ξ_s, d_s, k_s; 
        offset_size=0.15,
        φ=sphere_embedding,
        offsets=O,
        Ω=Ω
        )



# --------------------------------- simulate --------------------------------- #
dt = 0.5
duration = 500  # ms

nframes = (Int ∘ round)(duration / dt)
trajectory = Trajectory(can; T = nframes, 
        σθ=0.0, θ₀=deg2rad(90), σv=0.0, μv=.1 
)
simulation = Simulation(can, trajectory; b₀=1.0, η=0.0)

h = run_simulation(
    simulation,
    frame_every_n = 20,
    discard_first_ms = 0,
    average_over_ms = 20,
)
