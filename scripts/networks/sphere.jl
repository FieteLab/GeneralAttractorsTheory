using Distances
using Term

import GeneralAttractors: SphericalDistance
using GeneralAttractors
using GeneralAttractors.Kernels
using GeneralAttractors: lerp
using GeneralAttractors.ManifoldUtils
import GeneralAttractors.ManifoldUtils: sphere_embedding, ψx, ψy, ψz, fibonacci_sphere
import GeneralAttractors.Can: OneForm
import Manifolds: Sphere as 𝕊
import Manifolds: uniform_distribution

println(Panel("Creating sphere attractor", style="green", justify=:center))



# number of neurons
m =40
n = m^2

# get neurons on S² ⊂ ℝ³
sphere = 𝕊(2)
# X = hcat(rand(uniform_distribution(sphere), n)...)  # 3 × N neurons coordinate
X = fibonacci_sphere(n)

# get neurons indices
I = [(i, ) for i in 1:size(X, 2)]

# distance metric on the unit sphere
d_s = SphericalDistance()

# kernel  
k_s = LocalGlobalKernel(α = 0.5, σ = 0.2, β = 0.5)

# cover space
cover = CoverSpace(S²)  # trivial cover space

# define offset vector fields
offsets = [
    p -> ψx(p),
    p -> -ψx(p),
    p -> ψy(p),
    p -> -ψy(p),
    p -> ψz(p),
    p -> -ψz(p),
]
offset_size = .15

# define one forms
α = 1/offset_size .* 2
Ω = [
    OneForm(1, (x, y, z) -> α * ψx(x,y,z)),
    OneForm(2, (x, y, z) -> -α * ψx(x,y,z)),
    OneForm(3, (x, y, z) -> α * ψy(x,y,z)),
    OneForm(4, (x, y, z) -> -α * ψy(x,y,z)),
    OneForm(5, (x, y, z) -> α * ψz(x,y,z)),
    OneForm(6, (x, y, z) -> -α * ψz(x,y,z)),
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
    Ω=Ω
)
