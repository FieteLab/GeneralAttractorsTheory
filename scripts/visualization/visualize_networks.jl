using GeneralAttractors
using GeneralAttractors.Kernels
using GeneralAttractors.Manifolds
using GeneralAttractors: lerp
using Distances

import GeneralAttractors.Manifolds: sphere_embedding

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
    [cos(2π), sin(2π)],
]

Ω = OneForm[
    OneForm(1, x -> 1),
    OneForm(1, x -> 1),
    OneForm(1, x -> -1),
    OneForm(1, x -> -1),
    OneForm(2, x -> 1),
    OneForm(2, x -> 1),
    OneForm(2, x -> -1),
    OneForm(2, x -> -1),
]


sph = CAN("sphere", cover, n, ξ_s, d_s, k_s; 
        offset_size=0.15,
        φ=sphere_embedding,
        offsets=O,
        # Ω=Ω
        )




# show_connectivity(ring_attractor; plot_title = "Ring attractor connectivity")

# show_connectivity(torus_attractor; plot_title = "Torus attractor connectivity")

# show_connectivity(mobius_attractor; plot_title = "Mobius attractor connectivity")

show_connectivity(
    sph;
    xlabel = "longitude",
    ylabel = "latitude",
    plot_title = "Sphere attractor connectivity",
    idxs = [1, 30, 2100, 1500, 800],
    aspect_ratio=0.5,
    size=(800, 400)
)
