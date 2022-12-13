using Distances
using Term


using GeneralAttractors
using GeneralAttractors.Kernels
using GeneralAttractors: lerp
using GeneralAttractors.ManifoldUtils
import GeneralAttractors.ManifoldUtils: Ring, ring_ψ


println(Panel("Creating ring attractor", style = "green", justify = :center))

# neurons position and distance function
n = (200,)  # number of neurons in the ring

# neurons coordinates and metric
ξ_r(i::Int)::Vector = [lerp(i, n[1], 0.0, 2π - 2π / n[1])]  # neurons coordinates function
d_r = PeriodicEuclidean([2π])  # distance function

# kernel
k_r = LocalGlobalKernel(α = 2.5, σ = 0.25, β = 2.5)

# cover map
cover = CoverSpace(Ring())

# offsets and one forms
offset_size = .1
offsets = [
    p -> ring_ψ(p),
    p -> -ring_ψ(p)
]

Ω = OneForm[
    OneForm(1, (x) -> ring_ψ(x)),
    OneForm(1, (x) -> -ring_ψ(x))
]

# make network
ringcan = CAN("ring", cover, n, ξ_r, d_r, k_r; 
    offsets = offsets,
    Ω = Ω,
    offset_size = offset_size, 
    σ = :softrelu, 
    α = 50) # 120
