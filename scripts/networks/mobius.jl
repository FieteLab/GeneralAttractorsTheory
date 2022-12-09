using Distances
using Term


using GeneralAttractors
using GeneralAttractors.Kernels
using GeneralAttractors: lerp
using GeneralAttractors.ManifoldUtils
import GeneralAttractors.ManifoldUtils: Manifoldℝ², Mobius, ψ_t, ψ_θ1, ψ_θ2
import GeneralAttractors: MobiusEuclidean, mobius_embedding

println(Panel("Creating Mobius attractor", style = "green", justify = :center))

# number of neurons
n = ((Int ∘ round)(1 / 0.025), (Int ∘ round)(2π / 0.1))

# cover space
mfld = Mobius()
cover = CoverSpace(mfld)

# coordinates function (from neurons index to lattice coordintes)
ξ_m(i::Int, j::Int)::Vector =
    [lerp(i, n[1], -1 / 2, 1 / 2), lerp(j, n[2], 0.0, 2π - 2π / n[2])]  # ∈ [-1/2, 1/2] × [0, 2π]

# metric
d_m = MobiusEuclidean()

# connectivity kernel
# k_m = DiffOfExpKernel(; λ = 1.5)
k_m = LocalGlobalKernel(α = 0.25, σ = 0.1, β = 0.25)


# define offset vector fields
offsets =
    [p -> ψ_t(p), p -> -ψ_t(p), p -> ψ_θ1(p), p -> -ψ_θ1(p), p -> ψ_θ2(p), p -> -ψ_θ2(p)]
offset_size = 0.1

# define one forms
α = 1 / offset_size .* 2
Ω = [
    OneForm(1, (t, θ) -> α * ψ_t(t, θ)),
    OneForm(2, (t, θ) -> -α * ψ_t(t, θ)),
    OneForm(3, (t, θ) -> α * ψ_θ1(t, θ)),
    OneForm(4, (t, θ) -> -α * ψ_θ1(t, θ)),
    OneForm(5, (t, θ) -> α * ψ_θ2(t, θ)),
    OneForm(6, (t, θ) -> -α * ψ_θ2(t, θ)),
]



# construct network
mobiuscan = CAN(
    "mobius",
    cover,
    n,
    ξ_m,
    d_m,
    k_m;
    offset_size = offset_size,
    offsets = offsets,
    Ω = Ω,
)
