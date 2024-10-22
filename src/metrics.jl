using Distances: euclidean, Metric
import Distances
import LinearAlgebra: ⋅
import Manifolds: Sphere as 𝕊
import Manifolds: distance as mdist



"""
Sperical distance on the unit sphere by Manifolds.jl
"""
struct SphericalDistance <: Metric
    s::𝕊
end

SphericalDistance() = SphericalDistance(𝕊(2))

(dist::SphericalDistance)(p1, p2) = mdist(dist.s, p1, p2)
Distances.eval_op(::SphericalDistance, ::Float64, ::Float64) = 1.0




# struct MobiusEuclidean <: UnionMetric end

# """
# Compute the euclidean distnce between points in ℝ³
# """
# (dist::MobiusEuclidean)(p, q) = euclidean(
#         mobius_embedding(p), mobius_embedding(q)
#     )



"""
    MobiusEuclidean{W}

Euclidean metric on a Mobius strip space:

    ┏━━━━━━ << ━━━━━━┓
    ┃                ┃
    ┃                ┃
    ┃                ┃
    ┃                ┃
    ┃                ┃
    ┃                ┃
    ┃    ⨀           ┃
    ┗━━━━━━ >> ━━━━━━┛

It assumes a MB parametrized by:
    - t ∈ [-1/2, 1/2]
    - θ ∈ [0, 2π]

To compute the distance:
    1. get Δθ
    2. if Δθ < 2π/2 p,q are considered on the "same side" and
        their distance is just given by a PeriodicEuclidean metric
    3. if Δθ > 2π/2 we change p=(t, θ) to be p̂=(1-t, θ) and then
        use the PeriodicEuclidean metric

"""
struct MobiusEuclidean <: UnionMetric
    T::Float64  # "height" of the mfld in the non-periodic direction
    th::Float64  # treshold distance for points "on the same side"
    periodic::PeriodicEuclidean
end

MobiusEuclidean() = MobiusEuclidean(1.0, π, PeriodicEuclidean([Inf, 2π]))

function (dist::MobiusEuclidean)(q, p)
    @inbounds begin
        Δθ = abs(q[2] - p[2])
        ŷ = Δθ > dist.th ? [-p[1], p[2]] : p
        return dist.periodic(q, ŷ)
    end
end

# # Distances.result_type(::MobiusEuclidean, ::Float64, ::Float64) = Float64
Distances.eval_op(::MobiusEuclidean, ::Float64, ::Float64) = 1.0



# ---------------------------------------------------------------------------- #
#                                 klein bottle                                 #
# ---------------------------------------------------------------------------- #

struct KleinBottleEuclidean <: UnionMetric
    th::Float64  # threshold distance for points "on the same side"
    periodic::PeriodicEuclidean
end

KleinBottleEuclidean() = KleinBottleEuclidean(π, PeriodicEuclidean([2π, 2π]))

function (dist::KleinBottleEuclidean)(q, p)
    @inbounds begin
        Δx = abs(q[1] - p[1])
        Δy = abs(q[2] - p[2])
        
        # Handle x-direction (torus-like behavior)
        x̂ = Δx > π ? (q[1] > p[1] ? q[1] - 2π : q[1] + 2π) : q[1]
        
        # Handle y-direction (Möbius strip-like behavior)
        if Δy > dist.th
            ŷ = q[2] > p[2] ? q[2] - 2π : q[2] + 2π
            x̂ = 2π - x̂  # Flip x-coordinate when crossing y-boundary
        else
            ŷ = q[2]
        end
        
        q̂ = [x̂, ŷ]
        return dist.periodic(q̂, p)
    end
end

Distances.eval_op(::KleinBottleEuclidean, ::Float64, ::Float64) = 1.0
