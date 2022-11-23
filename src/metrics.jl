using Distances: UnionMetric, PeriodicEuclidean, euclidean, Metric
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

# (dist::SphericalDistance)(p1, p2) = mdist(dist.s, p1, p2)
(dist::SphericalDistance)(p1, p2) = begin
    c = max(p1 ⋅ p2, -1)
    c = min(c, 1)
    acos(c)
end
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

# """ Parametrized line in ℝ² """
# line(p, q, t) = [
#     t*p[1]+(1-t)*q[1],
#     t*p[2]+(1-t)*q[2],
# ]

# """ Length on a line in ℝ² """
# # ℓ(l::Matrix) = sum(
# #     sqrt.(diff(l[1, :]).^2 + diff(l[2, :]).^2)
# # )
# ℓ(l::Matrix) = sum(
#     sqrt.(diff(l[1, :]).^2 + diff(l[2, :]).^2 + diff(l[3, :]).^2)
# )

# function (dist::MobiusEuclidean)(p, q)
#     # first get a direct line between p,q
#     l1 = hcat(
#         map(
#             t->line(p, q, t), 0:.01:1
#         )...
#     )

#     # then get two lines segments going around the "other side"
#     l2_a = hcat(
#         map(
#             t->line(p, [0, q[2]-2π], t), 0:.01:1
#         )...
#     )

#     stopper = findfirst(l2_a[2, :] .>= 0.0)
#     l2_a = l2_a[:, stopper:end]
#     l2_b = hcat(
#         map(
#             t->line([-l2_a[1, 1], 2π], q, t), 0:.01:1
#         )...
#     )

#     # then get the smallest length between the two lines
#     return min((ℓ ∘ mobius_embedding)(l1), (ℓ ∘ mobius_embedding)(l2_a) + (ℓ ∘ mobius_embedding)(l2_b))
# end


function (dist::MobiusEuclidean)(q, p)
    @inbounds begin
        Δθ = abs(q[2] - p[2])
        ŷ = Δθ > dist.th ? [-p[1], p[2]] : p
        return dist.periodic(q, ŷ)
    end
end

# # Distances.result_type(::MobiusEuclidean, ::Float64, ::Float64) = Float64
Distances.eval_op(::MobiusEuclidean, ::Float64, ::Float64) = 1.0
