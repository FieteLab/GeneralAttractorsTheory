using Distances: euclidean, Metric
import Distances
import LinearAlgebra: ⋅
import Manifolds: Sphere as 𝕊
import Manifolds: distance as mdist
using Graphs, SimpleWeightedGraphs



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

using LinearAlgebra



struct KleinBottleEuclidean <: UnionMetric
    n_u::Int
    n_v::Int
    neighbors::Int
    graph::SimpleWeightedGraph
    distance_matrix::Matrix{Float64}
end


function klein_bottle_embedding(u, v)
    # Parametric equations for the 4D Klein bottle
    x1 = (2 + cos(v)) * cos(u)
    x2 = (2 + cos(v)) * sin(u)
    x3 = sin(v)
    x4 = sin(v) * cos(u / 2)
    return [x1, x2, x3, x4]  # 4D point
end


function create_klein_bottle_graph(n_u, n_v, neighbors)
    println("creating klein bottle graph")
    n_vertices = n_u * n_v
    g = SimpleWeightedGraph(n_vertices)
    
    # Embed all points in the embedding space
    points = [klein_bottle_embedding(2π * (i-1) / n_u, 2π * (j-1) / n_v) for i in 1:n_u, j in 1:n_v]
    
    for i in 1:n_u
        for j in 1:n_v
            vertex1 = (i-1) * n_v + j
            p1 = points[i, j]
            
            # Calculate distances to all other points
            distances = [(k, l, norm(p1 - points[k, l])) for k in 1:n_u, l in 1:n_v if (k, l) != (i, j)]
            
            # Sort by distance and take the 8 closest points
            closest_points = sort(distances, by=x->x[3])[1:neighbors]
            
            for (k, l, dist) in closest_points
                vertex2 = (k-1) * n_v + l
                if !has_edge(g, vertex1, vertex2)
                    add_edge!(g, vertex1, vertex2, dist)
                end
            end
        end
    end
    
    return g
end

function KleinBottleEuclidean(resolution=50, neighbors=10) 
    n_u, n_v = resolution, resolution
    graph = create_klein_bottle_graph(n_u, n_v, neighbors)
    distance_matrix = floyd_warshall_shortest_paths(graph).dists
    return KleinBottleEuclidean(n_u, n_v, neighbors, graph, distance_matrix)
end


function geodesic_distance_on_kb_graph(u1, v1, u2, v2, distance_matrix, n_u, n_v)
    start_vertex = clamp(round(Int, u1 / (2π) * n_u) * n_v + round(Int, v1 / (2π) * n_v) + 1, 1, n_u * n_v)
    end_vertex = clamp(round(Int, u2 / (2π) * n_u) * n_v + round(Int, v2 / (2π) * n_v) + 1, 1, n_u * n_v)
    
    return distance_matrix[start_vertex, end_vertex]
end

function (dist::KleinBottleEuclidean)(q, p)
    u1, v1 = q
    u2, v2 = p
    return geodesic_distance_on_kb_graph(u1, v1, u2, v2, dist.distance_matrix, dist.n_u, dist.n_v)
end

Distances.eval_op(::KleinBottleEuclidean, ::Float64, ::Float64) = 1.0
