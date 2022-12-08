using Distances
using Term


using GeneralAttractors
using GeneralAttractors.Kernels
using GeneralAttractors: lerp
using GeneralAttractors.ManifoldUtils
import GeneralAttractors.ManifoldUtils: Ring


println(Panel("Creating ring attractor", style = "green", justify = :center))

# neurons position and distance function
n = (256,)  # number of neurons in the ring

# neurons coordinates and metric
ξ_r(i::Int)::Vector = [lerp(i, n[1], 0.0, 2π-2π/n[1])]  # neurons coordinates function
d_r = PeriodicEuclidean([2π])  # distance function

# kernel
# k_r = DiffOfExpKernel(; λ = 0.1, β=1.5)
k_r = LocalGlobalKernel(α = 0.5, σ = 5.0, β = 0.5)

# cover map
cover = CoverSpace(Ring())

# make network
ringcan = CAN(
    "ring",
    cover,
    n,
    ξ_r,
    d_r,
    k_r;
    offset_size = 1.0,
    σ = :softrelu,
    α= 0.9
)
