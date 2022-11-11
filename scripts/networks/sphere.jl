using Distances
using Term


using GeneralAttractors
using GeneralAttractors.Kernels
using GeneralAttractors: lerp
using GeneralAttractors.ManifoldUtils
import GeneralAttractors.ManifoldUtils: sphere_embedding, ψx, ψy, ψz, ψxS², ψyS², ψzS²
import GeneralAttractors.Can: OneForm

println(Panel("Creating sphere attractor", style="green", justify=:center))


# number of neurons
n = (64, 64)

""" coordinates map on the sphere """
function ξ_s(i::Int, j::Int)::Vector
    [lerp(i, n[1], -π, π), lerp(j, n[2], -π / 2, π / 2)]
end

# distance metric
d_s = SphericalAngle()

# kernel
# k_s = DiffOfExpKernel(; λ = 0.75)
k_s = MexicanHatKernel(α=.015, σ=.3, β=0.01)

# cover space
cover = CoverSpace(S², S², (x, y) -> [x, y])

# define offset vector fields
offsets = [
    ψx,
    p -> -ψx(p),
    ψy,
    p -> -ψy(p),
    ψz,
    p -> -ψz(p)]

# define one forms
α=1
Ω = [
    OneForm(1, (x, y) -> α * ψxS²([x,y])),
    OneForm(2, (x, y) -> -α * ψxS²([x,y])),
    OneForm(3, (x, y) -> α * ψyS²([x,y])),
    OneForm(4, (x, y) -> -α * ψyS²([x,y])),
    OneForm(5, (x, y) -> α * ψzS²([x,y])),
    OneForm(6, (x, y) -> -α * ψzS²([x,y])),
]


# construct CAN
spherecan = CAN(
    "sphere",
    cover,
    n,
    ξ_s,
    d_s,
    k_s;
    offset_size = 0.25,
    offsets = offsets,
    φ=sphere_embedding,
    Ω=Ω
)
