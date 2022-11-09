using GeneralAttractors
using GeneralAttractors.Kernels
using GeneralAttractors.Manifolds
using GeneralAttractors: lerp
using Distances


n = (64, 64)
function ξ_s(i::Int, j::Int)::Vector
    [lerp(i, n[1], -π, π), lerp(j, n[2], -π / 2, π / 2)]
end
d_s = SphericalAngle()
k_s = DiffOfExpKernel(; λ = 0.75)

cover = CoverSpace(S², S², (x, y) -> [mod(x, 1), mod(y, 1)])
sph = CAN("sphere", cover, n, ξ_s, d_s, k_s; offset_size=[0.2, 0.2, 0.1, 0.1])




# show_connectivity(ring_attractor; plot_title = "Ring attractor connectivity")

# show_connectivity(torus_attractor; plot_title = "Torus attractor connectivity")

# show_connectivity(mobius_attractor; plot_title = "Mobius attractor connectivity")

show_connectivity(
    sph;
    xlabel = "longitude",
    ylabel = "latitude",
    plot_title = "Sphere attractor connectivity",
)
