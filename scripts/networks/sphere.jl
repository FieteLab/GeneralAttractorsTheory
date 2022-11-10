using Distances

using GeneralAttractors
using GeneralAttractors.Kernels
using GeneralAttractors: lerp
using GeneralAttractors.ManifoldUtils
import GeneralAttractors.ManifoldUtils: sphere_embedding

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


"""
weight offsets vector field. 
These are defined as rotational vector fields
on the unit sphere in ℝ³ with rotations around X, Y, Z.
These are killing fields and are divergence free.


Given a point p=(x, y, z) ∈ S² ⊂ ℝ²,
get a vector (fx(p)∂x, fy(p)∂y, fz(p)∂z)
tangent to the sphere and correpsonding to a rotation.
"""

∂x = [1, 0, 0]
∂y = [0, 1, 0]
∂z = [0, 0, 1]

""" rotation around X axis """
ψx(x, y, z) = z * ∂y - y * ∂z
ψx(p) = ψx(p...) .- p

""" rotation around Y axis """
ψy(x, y, z) = z * ∂x - x * ∂z
ψy(p) = ψy(p...) .- p

""" rotation around Z axis """
ψz(x, y, z) = x * ∂y - y * ∂x
ψz(p) = ψz(p...) .- p

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
