using Distances
using Term

import GeneralAttractors: SphericalDistance
using GeneralAttractors
using GeneralAttractors.Kernels
using GeneralAttractors: lerp
using GeneralAttractors.ManifoldUtils
import GeneralAttractors.ManifoldUtils: sphere_embedding, ψx, ψy, ψz, ψxS², ψyS², ψzS²
import GeneralAttractors.Can: OneForm
import Manifolds: Sphere as 𝕊
import Manifolds: uniform_distribution

println(Panel("Creating sphere attractor", style="green", justify=:center))


"""
(almost) equally spaced points on the sphere.
https://stackoverflow.com/questions/9600801/evenly-distributing-n-points-on-a-sphere
"""
function fibonacci_sphere(n=1000)
    points = zeros(3 ,n)
    ϕ = π * (3 - √5)  # golden angle in radians

    for i in 1:n
        y = 1 - (i / float(n - 1)) * 2  # y goes from 1 to -1
        radius = √(Complex(1 - y * y)) |> real  # radius at y

        θ = ϕ * i  # golden angle increment

        x = cos(θ) * radius
        z = sin(θ) * radius

        points[:, i] = [x, y, z]
    end
    return points
end


# number of neurons
m = 32
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
# k_s = DiffOfExpKernel(; λ = 5.0) 
k_s = LocalGlobalKernel(α = 0.5, σ = 1.0, β = 0.6)

# cover space
cover = CoverSpace(S², S², (x, y) -> [x, y])

# define offset vector fields
offsets = [
    p -> ψx(p),
    p -> -ψx(p),
    p -> ψy(p),
    p -> -ψy(p),
    p -> ψz(p),
    p -> -ψz(p),
]
offset_size = .1

# define one forms
α = 1/offset_size .* 100
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
