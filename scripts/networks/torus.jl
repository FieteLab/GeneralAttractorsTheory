using Distances
using Term


using GeneralAttractors
using GeneralAttractors.Kernels
using GeneralAttractors: lerp
using GeneralAttractors.ManifoldUtils


println(Panel("Creating torus attractor", style="green", justify=:center))


n = (32, 32)

cover = CoverSpace(ℝ², T, (x, y) -> [mod(x - 32, 64), mod(y - 32, 64)])
function ξ_t(i::Int, j::Int)::Vector  # neurons coordinates function
    n̂_i, n̂_j = Int(n[1] / 2), Int(n[2] / 2)
    [lerp(i, n[1], -n̂_i, n̂_i), lerp(j, n[2], -n̂_j, n̂_j)]   # ∈ [-n/2, n/2] × [-n/2, n/2]
end
d_t = PeriodicEuclidean([n...])  # distance function over a torus manifold

# connectivity kernel 
# k_t = DiffOfExpKernel(; λ = 13.0)
k_t = LocalGlobalKernel(α = 0.25, σ = 25.0, β = 0.25)

# one forms
Ω = OneForm[
    OneForm(1, x -> sin(2x / n[1]) + 1.25),
    OneForm(1, x -> -(sin(2x / n[1]) + 1.25)),
    OneForm(2, x -> sin(2x / n[2]) + 1.25),
    OneForm(2, x -> -(sin(2x / n[2]) + 1.25)),
]

# make network
toruscan = CAN("torus", cover, n, ξ_t, d_t, k_t; 
    offset_size=1.0, 
    # Ω = Ω
    )
