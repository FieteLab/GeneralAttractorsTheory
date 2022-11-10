using Distances

using GeneralAttractors
using GeneralAttractors.Kernels
using GeneralAttractors: lerp
using GeneralAttractors.ManifoldUtils
import GeneralAttractors.ManifoldUtils: sphere_embedding, ψx, ψy, ψz

# number of neurons
n = (64, 64)

""" coordinates map on the sphere """
function ξ_s(i::Int, j::Int)::Vector
    [lerp(i, n[1], -π, π), lerp(j, n[2], -π / 2, π / 2)]
end

# distance metric
d_s = SphericalAngle()

# kernel
k_s = DiffOfExpKernel(; λ = 0.75)

# cover space
cover = CoverSpace(S², S², (x, y) -> [mod(x, 1), mod(y, 1)])



# construct CAN
spherecan = CAN(
    "sphere",
    cover,
    n,
    ξ_s,
    d_s,
    k_s;
    offset_size = 0.25,
    offsets = [ψx, p -> -ψx(p), ψy, p -> -ψy(p), ψz, p -> -ψz(p)],
    φ=sphere_embedding
)
