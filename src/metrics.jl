
using Distances: UnionMetric, PeriodicEuclidean, euclidean
import Distances

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
    - if the distance |Δx|>π (in the periodic direction)
        - define p̂ = (x₂, 1-y₂)   and compute d(q, p̂)
            where d = PeriodicEuclidean([2π, Inf])
        - otherwise the distance is d(q, p)

Δx is computed using a PeriodicEuclidean([2π]) metric. 

"""
struct MobiusEuclidean <: UnionMetric
    periodic1D::PeriodicEuclidean
    periodic::PeriodicEuclidean
end

MobiusEuclidean() = MobiusEuclidean(
        PeriodicEuclidean([2π]),
        PeriodicEuclidean([2π, Inf])
    )


function (dist::MobiusEuclidean)(x, y)
    @inbounds begin
        Δx = abs((x[1] - y[1]))

        if Δx > π
            return dist.periodic(
                x, [y[1], 1-y[2]]
            )
        else
            return dist.periodic(x, y)
        end
    end
end

# Distances.result_type(::MobiusEuclidean, ::Float64, ::Float64) = Float64
Distances.eval_op(::MobiusEuclidean, ::Float64, ::Float64) = 1.0

