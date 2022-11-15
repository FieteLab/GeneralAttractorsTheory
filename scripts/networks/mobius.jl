using Distances
using Term


using GeneralAttractors
using GeneralAttractors.Kernels
using GeneralAttractors: lerp
using GeneralAttractors.ManifoldUtils
import GeneralAttractors.ManifoldUtils: Manifoldℝ², Mobius
import GeneralAttractors: MobiusEuclidean

println(Panel("Creating Mobius attractor", style="green", justify=:center))

# number of neurons
m = 32
n = (m, m)

# cover space
mfld = Mobius(2π)
cover = CoverSpace(mfld)

# coordinates function (from neurons index to lattice coordintes)
ξ_m(i::Int, j::Int)::Vector = [
    lerp(i, n[1], 0.0, 2π), 
    lerp(j, n[2], 0.0, 2π)
]  # ∈ [0, 2π] × [0, 2π]

# metric
d_m = MobiusEuclidean(2π)

# connectivity kernel
# k_m = DiffOfExpKernel(; λ = 1.5)
k_m = LocalGlobalKernel(α = 0.25, σ = 1.0, β = 0.25)

# construct network
mobius_attractor = CAN("mobius", 
        cover, n, ξ_m, d_m, k_m;
        offset_size=0.2,
)
