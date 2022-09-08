
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
    - define p̂ = (x₂, 1-y₂)  if the shortest path p̄q̄ is
        through the "glued" part of the strip. Otherwise
        p̂ = p.
    - compute the distance with a PeriodicEuclidean metric


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
    # flip y-coord for odd trips around the strip
    @inbounds begin
        Δx = dist.periodic1D(x[1], y[1])
        if Δx < 2π - Δx
            return dist.periodic(x, y)
        else
            ŷ = [y[1], 1-y[2]]
        end
    end

    # get distance
    return dist.periodic(x, ŷ)
end

# Distances.result_type(::MobiusEuclidean, ::Float64, ::Float64) = Float64
Distances.eval_op(::MobiusEuclidean, ::Float64, ::Float64) = 1.0

