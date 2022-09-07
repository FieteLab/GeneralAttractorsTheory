
using Distances: UnionMetric, PeriodicEuclidean, euclidean


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
    periodic::PeriodicEuclidean
end

MobiusEuclidean() = MobiusEuclidean(PeriodicEuclidean([2π, Inf]))


function (dist::MobiusEuclidean)(x, y)
    @assert length(x) == length(y) == 2 "Inputs to Mobius Distance should be 2-vectors"


    # flip y-coord for odd trips around the strip
    if abs(x[1]-y[1]) < abs(x[1]-(2π-y[1]))
        ŷ = y
    else
        @inbounds ŷ = [y[1], 1-y[2]]
    end
    # @inbounds ŷ = [y[1], 1-y[2]]

    # get distance
    return dist.periodic(x, ŷ)
end
