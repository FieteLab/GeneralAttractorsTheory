using Distances
using Term

import GeneralAttractors: SphericalDistance
using GeneralAttractors
using GeneralAttractors.Kernels
using GeneralAttractors: lerp
using GeneralAttractors.ManifoldUtils
import GeneralAttractors.ManifoldUtils: sphere_embedding, ψx, ψy, ψz, fibonacci_sphere
import GeneralAttractors.Can: OneForm


println(Panel("Creating sphere attractor", style = "green", justify = :center))



# number of neurons
m = 40
n = m^2

# get neurons on S² ⊂ ℝ³
X = fibonacci_sphere(n)

# get neurons indices
I = [(i,) for i = 1:size(X, 2)]

# distance metric on the unit sphere
d_s = SphericalDistance()
# d_s = SphericalAngle()

# kernel  
k_s = LocalGlobalKernel(α = 2.5, σ = 40.5, β = 2.5)

# cover space
cover = CoverSpace(S²)  # trivial cover space

# define offset vector fields
offsets = [p -> ψx(p), p -> -ψx(p), p -> ψy(p), p -> -ψy(p), p -> ψz(p), p -> -ψz(p)]
offset_size = 0.1

# define one forms
Ω = [
    OneForm(1, (x, y, z) -> offset_size * ψx(x, y, z)),
    OneForm(1, (x, y, z) -> -offset_size * ψx(x, y, z)),
    OneForm(2, (x, y, z) -> offset_size * ψy(x, y, z)),
    OneForm(2, (x, y, z) -> -offset_size * ψy(x, y, z)),
    OneForm(3, (x, y, z) -> offset_size * ψz(x, y, z)),
    OneForm(3, (x, y, z) -> -offset_size * ψz(x, y, z)),
]


# construct CAN
spherecan = CAN(
    "sphere",
    cover,
    (m, m),
    I,
    X,
    d_s,
    k_s;
    offset_size = offset_size,
    offsets = offsets,
    Ω = Ω,
    α = 46,
)
