using LinearAlgebra, DifferentialEquations, ForwardDiff, Plots
using Graphs, SimpleWeightedGraphs

function klein_bottle_embedding(u, v)
    # Parametric equations for the 4D Klein bottle
    x1 = (2 + cos(v)) * cos(u)
    x2 = (2 + cos(v)) * sin(u)
    x3 = sin(v)
    x4 = sin(v) * cos(u / 2)
    return [x1, x2, x3, x4]  # 4D point
end

# function klein_bottle_embedding(u, v)
#     # Parametric equations for a 3D self-intersecting Klein bottle
#     x = (2 + cos(u / 2) * sin(v) - sin(u / 2) * sin(2 * v)) * cos(u)
#     y = (2 + cos(u / 2) * sin(v) - sin(u / 2) * sin(2 * v)) * sin(u)
#     z = sin(u / 2) * sin(v) + cos(u / 2) * sin(2 * v)
#     return [x, y, z]  # 3D point
# end



function create_klein_bottle_graph(n_u, n_v)
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
            closest_points = sort(distances, by=x->x[3])[1:8]
            
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

function geodesic_distance(u1, v1, u2, v2, graph, n_u, n_v)
    start_vertex = round(Int, u1 / (2π) * n_u) * n_v + round(Int, v1 / (2π) * n_v) + 1
    end_vertex = round(Int, u2 / (2π) * n_u) * n_v + round(Int, v2 / (2π) * n_v) + 1
    
    path = dijkstra_shortest_paths(graph, start_vertex)
    return path.dists[end_vertex]
end

function plot_klein_bottle(n_u=50, n_v=50)
    u = range(0, 2π, length=n_u)
    v = range(0, 2π, length=n_v)
    
    points = [klein_bottle_embedding(u, v) for u in u, v in v]
    x = [p[1] for p in points]
    y = [p[2] for p in points]
    z = [p[3] for p in points]
    
    p = scatter3d(x, y, z, markersize=1, color=:viridis, alpha=0.7, legend=false, label=nothing)
    # set axes limits
    xlims!(p, -3, 3)
    ylims!(p, -3, 3)
    zlims!(p, -3, 3)
end

function plot_graph_on_klein_bottle(graph, n_u, n_v)
    p = plot()
    
    for edge in edges(graph)
        source = edge.src
        dest = edge.dst
        u1, v1 = 2π * ((source-1) ÷ n_v) / n_u, 2π * ((source-1) % n_v) / n_v
        u2, v2 = 2π * ((dest-1) ÷ n_v) / n_u, 2π * ((dest-1) % n_v) / n_v
        p1 = klein_bottle_embedding(u1, v1)
        p2 = klein_bottle_embedding(u2, v2)
        plot!(p, [p1[1], p2[1]], [p1[2], p2[2]], [p1[3], p2[3]], color=:black, linewidth=1, alpha=0.3, label=nothing)
    end
    # set axes limits
    xlims!(p, -3, 3)
    ylims!(p, -3, 3)
    zlims!(p, -3, 3)
    p
end

function plot_geodesic(u1, v1, u2, v2, graph, n_u, n_v)
    start_vertex = round(Int, u1 / (2π) * n_u) * n_v + round(Int, v1 / (2π) * n_v) + 1
    end_vertex = round(Int, u2 / (2π) * n_u) * n_v + round(Int, v2 / (2π) * n_v) + 1
    
    path = dijkstra_shortest_paths(graph, start_vertex)
    shortest_path = enumerate_paths(path, end_vertex)
    
    # p = plot_klein_bottle()
    p = plot_graph_on_klein_bottle(graph, n_u, n_v)
    
    for i in 1:length(shortest_path)-1
        source, dest = shortest_path[i], shortest_path[i+1]
        u1, v1 = 2π * ((source-1) ÷ n_v) / n_u, 2π * ((source-1) % n_v) / n_v
        u2, v2 = 2π * ((dest-1) ÷ n_v) / n_u, 2π * ((dest-1) % n_v) / n_v
        p1 = klein_bottle_embedding(u1, v1)
        p2 = klein_bottle_embedding(u2, v2)
        plot!(p, [p1[1], p2[1]], [p1[2], p2[2]], [p1[3], p2[3]], color=:red, linewidth=2, label=nothing)
    end

    # set axes limits
    xlims!(p, -3, 3)
    ylims!(p, -3, 3)
    zlims!(p, -3, 3)
    
    p
end

# Example usage
n_u, n_v = 50, 25  # Reduced for faster visualization
klein_bottle_graph = create_klein_bottle_graph(n_u, n_v)

# Plot the Klein bottle
# p1 = plot_klein_bottle()
title!(p1, "Klein Bottle Surface")

# Plot the graph on the Klein bottle
p2 = plot_graph_on_klein_bottle(klein_bottle_graph, n_u, n_v)
title!(p2, "Graph on Klein Bottle")

# Compute and plot geodesic
u1, v1 = π, π-0.1
u2, v2 = 0, 0
distance = geodesic_distance(u1, v1, u2, v2, klein_bottle_graph, n_u, n_v)
p3 = plot_geodesic(u1, v1, u2, v2, klein_bottle_graph, n_u, n_v)
title!(p3, "Geodesic on Klein Bottle")
println("Distance: $distance")

# Display all plots
# plot(p1, p2, p3, layout=(1,3), size=(1200,400), legend=false, label=nothing)
display(plot(p2, p3))
