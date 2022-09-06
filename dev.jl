using GeneralAttractors

import Base.Iterators: product as ×  # cartesian product
import Distances: Metric, PeriodicEuclidean, pairwise
import StaticArrays: SVector, SA_F64, SMatrix
using Term.Progress
using Plots



using Logging
import Term.Logs: TermLogger
Logging.min_enabled_level(::TermLogger) = Logging.Debug






# ----------------------------------- ring ----------------------------------- #
@info "doing ring"
n = (100,)  # number of neurons in the ring
ξ_r(i::Int)::Vector = [2π/n[1]*i]  # neurons coordinates function
d = PeriodicEuclidean([2π])  # distance function


a = 1.0
α = 0.10315
l = 1.0
λ = 13
β = 3/(λ^2)
γ = 1.05 * β

k(x::Float64)::Float64 = a * exp(-γ*abs(x)^2) - exp(-β*abs(x)^2)

# @time can(n, ξ_r, d, Kernel(k));
# @time can(n, ξ_r, d, Kernel(k));

# ----------------------------------- torus ---------------------------------- #
@info "doing torus"
n = (64, 64)  # number of neurons in the ring
ξ_t(i::Int, j::Int)::Vector = begin # neurons coordinates function
    p_i, p_j = (i-1)/(n[1]-1), (j-1)/(n[2]-1) # ∈ [0, 1]
    [
        -π/2*(1-p_i)+π/2*(p_i), 
        -π/2*(1-p_j)+π/2*(p_j)
    ]  # ∈ [-π/2, π/2] × [-π/2, π/2]
end
d = PeriodicEuclidean([2π, 2π])  # distance function over a torus manifold

# can = @time CAN(n, ξ_t, d, Kernel(k))
show_connectivity(can,)