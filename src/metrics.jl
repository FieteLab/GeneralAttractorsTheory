using Distances: UnionMetric, PeriodicEuclidean, euclidean
import Distances
import LinearAlgebra: ⋅

import .ManifoldUtils: sphere_embedding

"""
Sperical distance on the unit sphere
https://mathworld.wolfram.com/SphericalDistance.html
"""
struct SphericalDistance <: UnionMetric end

function (dist::SphericalDistance)(x₁, x₂)
    @inbounds begin
        # get points on the embedded sphere manifold
        p₁ = sphere_embedding(x₁...)
        p₂ = sphere_embedding(x₂...)

        # get the distance
        return acos(p₁ ⋅ p₂)
    end
end

Distances.eval_op(::SphericalDistance, ::Float64, ::Float64) = 1.0



"""
    MobiusEuclidean{W}

Euclidean metric on a Mobius strip space:

    ┏━━━━━━━━━━━━━━━━━━┓
    x                  ┃
    ┃                  ┃
    ┃    ⨀             x
    ┗━━━━━━━━━━━━━━━━━━┛

The two `x` points are identified on the Mobius strip.
To compute the distance:
    - given two points q=(x₁, y₁), p=(x₂, y₂)
    - if the distance |Δx|>L/2 (half length in the periodic direction)
        - define p̂ = (x₂, 1-y₂)   and compute d(q, p̂)
            where d = PeriodicEuclidean([2π, Inf])
        - otherwise the distance is d(q, p)

Δx is computed using a PeriodicEuclidean([2π]) metric. 

"""
struct MobiusEuclidean <: UnionMetric
    x₀::Float64  # "height" of the mfld in the non-periodic direction
    periodic1D::PeriodicEuclidean
    periodic::PeriodicEuclidean
end

MobiusEuclidean() = MobiusEuclidean(1.0)

MobiusEuclidean(x₀::Float64) =
    MobiusEuclidean(x₀, PeriodicEuclidean([2π]), PeriodicEuclidean([2π, Inf]))


function (dist::MobiusEuclidean)(x, y)
    @inbounds begin
        Δx = abs(x[1] - y[1])
        th = dist.periodic1D.periods[1] / 2

        if Δx > th
            return dist.periodic(x, [y[1], dist.x₀ - y[2]])
        else
            return dist.periodic(x, y)
        end
    end
end

# Distances.result_type(::MobiusEuclidean, ::Float64, ::Float64) = Float64
Distances.eval_op(::MobiusEuclidean, ::Float64, ::Float64) = 1.0
