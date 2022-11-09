using GeneralAttractors
using GeneralAttractors.Kernels
using GeneralAttractors.Manifolds
using GeneralAttractors: lerp
using Distances

import GeneralAttractors.Manifolds: sphere_embedding
import GeneralAttractors: SphericalDistance

n = (64, 64)
function ξ_s(i::Int, j::Int)::Vector
    [lerp(i, n[1], -π, π), lerp(j, n[2], -π / 2, π / 2)]
end
d_s = SphericalDistance()
k_s = DiffOfExpKernel(; λ = 0.75)

cover = CoverSpace(S², S², (x, y) -> [x, y])
sph = CAN("sphere", cover, n, ξ_s, d_s, k_s; 
        offset_size=0.15,
        φ=sphere_embedding,
        )


show_connectivity(
    sph;
    xlabel = "longitude",
    ylabel = "latitude",
    plot_title = "Sphere attractor connectivity",
    idxs = [1, 30, 2100, 1500, 800],
    aspect_ratio=0.5,
    size=(800, 400)
)
