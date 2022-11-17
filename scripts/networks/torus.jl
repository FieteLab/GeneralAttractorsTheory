using Distances
using Term


using GeneralAttractors
using GeneralAttractors.Kernels
using GeneralAttractors: lerp
using GeneralAttractors.ManifoldUtils
import GeneralAttractors.ManifoldUtils: Manifoldℝ², Torus


println(Panel("Creating torus attractor", style="green", justify=:center))

# number of neurons
m = 40 # number of neurons in each dimension
n = (m, m) # number of neurons per dimension

# ℝ² → T cover map.
D = 300 # ℝ² considered to go from -D to D in each direction


""" \rho """
ρ(x, y) = [mod(x - m/2, m)-m/2, mod(y - m/2, m)-m/2] # ρ: ℝ² → T

"""
    ρⁱ(x, y)

Inverse of the cover map rho over the domain.
Given a point (x,y) in N it gives a set of points (x̂, ŷ)
in the cover space such that ρ(x̂, ŷ)=(x,y)
"""
function ρⁱ(x, y; n=20)
    pts = zeros(2, n^2)
    for (c, i) in enumerate(-n/2:(n/2-1)), (k, j) in enumerate(-n/2:(n/2-1))
        x̂ = x + m*i
        ŷ = y + m*j
        pts[:, (c-1)*n+k] = [x̂, ŷ]
    end
    return pts
end

cover = CoverSpace(
    Manifoldℝ²(D), 
    Torus(m), 
    ρ, ρⁱ
)

# define a function to get the coordinates of each neuron in the lattice
function ξ_t(i::Int, j::Int)::Vector  # neurons coordinates function
    n̂_i, n̂_j = Int(n[1] / 2), Int(n[2] / 2)
    [lerp(i, n[1], -n̂_i, n̂_i-1), lerp(j, n[2], -n̂_j, n̂_j-1)]   # ∈ [-n/2, n/2] × [-n/2, n/2]
    # [lerp(i, n[1], 0, 2π), lerp(j, n[2], 0, 2π)]   # ∈ [-n/2, n/2] × [-n/2, n/2]
end

# select a distance metric
d_t = PeriodicEuclidean([n...])  # distance function over a torus manifold
# d_t = PeriodicEuclidean([2π, 2π])  # distance function over a torus manifold

# connectivity kernel 
# k_t = DiffOfExpKernel(; λ = 13.0)
k_t = LocalGlobalKernel(α = 0.5, σ = 50.0, β = 0.5)

# one forms
# Ω = OneForm[
#     OneForm(1, x -> sin(2x / n[1]) + 1.25),
#     OneForm(1, x -> -(sin(2x / n[1]) + 1.25)),
#     OneForm(2, x -> sin(2x / n[2]) + 1.25),
#     OneForm(2, x -> -(sin(2x / n[2]) + 1.25)),
# ]

# make network
toruscan = CAN("torus", cover, n, ξ_t, d_t, k_t; 
offset_size=1.0,
    # Ω = Ω
    )
