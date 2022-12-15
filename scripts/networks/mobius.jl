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
n = ((Int ∘ round)(1 / 0.05), (Int ∘ round)(2π / 0.2))

# cover space
mfld = Mobius()
cover = CoverSpace(mfld)

# coordinates function (from neurons index to lattice coordintes)
ξ_m(i::Int, j::Int)::Vector =
    [
        lerp(i, n[1], mfld.xmin[1], mfld.xmax[1]), 
        lerp(j, n[2], mfld.xmin[2], mfld.xmax[2] - mfld.xmax[2] / n[2])
    ]  

# metric
d_m = MobiusEuclidean()

# connectivity kernel
k_m = LocalGlobalKernel(α = 2.5, σ = 1.5, β = 2.5)


# define offset vector fields
offsets =
    [
        p -> ψ_t(p),
        p -> -ψ_t(p),
        p -> ψ_θ1(p),
        p -> -ψ_θ1(p),
        p -> ψ_θ2(p), 
        p -> -ψ_θ2(p)
    ]
offset_size = 0.2

# define one forms
Ω = [
    OneForm(1, (t, θ) -> offset_size * ψ_t(t, θ)),
    OneForm(2, (t, θ) -> -offset_size * ψ_t(t, θ)),
    OneForm(3, (t, θ) -> offset_size * ψ_θ1(t, θ)),
    OneForm(4, (t, θ) -> -offset_size * ψ_θ1(t, θ)),
    OneForm(5, (t, θ) -> offset_size * ψ_θ2(t, θ)),
    OneForm(6, (t, θ) -> -offset_size * ψ_θ2(t, θ)),
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
    σ=:softrelu,
    α = 42,
)
