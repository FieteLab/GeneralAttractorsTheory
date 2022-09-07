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
# neurons position and distance function
n = (256,)  # number of neurons in the ring
ξ_r(i::Int)::Vector = [(i-1)/(n[1]-1)*2π]  # neurons coordinates function
d = PeriodicEuclidean([2π])  # distance function

# connectivity kernel
α, σ = .33, .25
mexicanhat(t) = α * 2/(√(3σ)*π^(1/4))*(1 - (t/σ)^2) * exp(-t^2/(2σ))
k(x::Float64)::Float64 = mexicanhat(x)

# construct network
can = CAN(n, ξ_r, d, Kernel(k); offset_strength=0.25);
show_connectivity(can) |> display

# ----------------------------------- torus ---------------------------------- #
@info "doing torus"

# neurons position and distance function
n = (48, 48)  # number of neurons in the ring

ξ_t(i::Int, j::Int)::Vector = begin # neurons coordinates function
    p_i, p_j = (i-1)/(n[1]-1), (j-1)/(n[2]-1) # ∈ [0, 1]
    [
        2π*(p_i), 
        2π*(p_j)
    ]  # ∈ [-π/2, π/2] × [-π/2, π/2]
end
d = PeriodicEuclidean([2π, 2π])  # distance function over a torus manifold

# connectivity kernel
k(x::Float64)::Float64 = begin    
    a = 1.0
    α = 0.10315
    l = 1.0
    λ = 13/2π
    β = 5/(λ^2)
    γ = 1.05 * β


    _x = abs(x)^2
    a * exp(-γ*_x) - exp(-β*_x)
end

# construct network
can = @time CAN(n, ξ_t, d, Kernel(k); offset_strength=.5)
show_connectivity(can,)

# TODO avoid having to bake kernel params in.


