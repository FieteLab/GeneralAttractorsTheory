using Plots
using Measures
using Distances: PeriodicEuclidean, evaluate, UnionMetric
import Base.Iterators: product as ×  # cartesian product


# ---------------------------------- kernel ---------------------------------- #
function Plots.plot(K::Kernel; kwargs...)
    x = -10:.1:10 |> collect
    y = K.(x)
    Plots.plot(
        x, 
        y,
        xlabel="distance Δx",
        ylabel="w",
        lw=2,
        label=nothing,
        color="black"; kwargs...)
end



# ---------------------------- distance functions ---------------------------- #
function plot_distance_1d(d::PeriodicEuclidean; kwargs...)
    x = 0:.1:d.periods[1] |> collect

    p = plot(; xlabel="x", ylabel="distance", kwargs...)
    colors = [:black, :red, :green, :blue]
    for (i, point) in enumerate(rand(x, 4))
        y = [evaluate(d, point, xx) for xx in x]
        plot!(x, y, lw=3, color=colors[i], label=nothing)
        vline!([point], color=colors[i], label=nothing)
    end
    display(p)
end


"""
    plot_distance_2d(
        d::UNion{PeriodicEuclidean, MobiusEuclidean}, 
        X::Vector{Vector})

Plot a 2D distance metric `d` given the set of 
points coordinates `x,y
"""
function plot_distance_2d(
    d, 
    x::Vector,
    y::Vector;
    points=nothing,
    kwargs...
    )

    # create lattice of points
    X = (x × y) |> collect
    X = [[x...] for x in vec(X)]

    # get distance from a point
    points = isnothing(points) ? rand(X, 4) : points
    pts = []
    for p in points
        Δx = [evaluate(d, p, x) for x in X]

        _plot=contourf(
            x, 
            y,
            reshape(Δx, length(x), length(y))',
            aspect_ratio=:equal,
            xlabel="x", ylabel="y",
            linewidth=0.25,
            grid=false,
            )
        scatter!([p[1]], [p[2]], color=:green, label=nothing, ms=10, alpha=1)
        push!(pts, _plot)
    end

    plt = plot(pts..., size=(600, 600); kwargs...)
    display(plt)
end

function plot_distance_2d(d::PeriodicEuclidean; kwargs...)
    upperbound(x, ) = isfinite(x) ? x : 1
    # get coordinates mesh
    x = range(0, upperbound(d.periods[1]), length=100) |> collect
    y = range(0, upperbound(d.periods[2]), length=100) |> collect

    plot_distance_2d(d, x, y; kwargs...)
end


function plot_distance_2d(d::MobiusEuclidean; kwargs...)
    x = 0:.075:2π |> collect
    y = 0:.075:1 |> collect
    X = (x × y) |> collect
    X = [[x...] for x in vec(X)]

    plot_distance_2d(d, x, y; points = [
        [0, 0],
        [3, 0],
        [.2, .5],
        [2π, 0]
    ], kwargs...)
end


"""
    plot_distance_function(d::PeriodicEuclidean)

Visualize N=dimensional periodic distance functions.
"""
function plot_distance_function(d::PeriodicEuclidean; kwargs...)
    ndims = length(d.periods)

    if ndims == 1
        plot_distance_1d(d; kwargs...)
    elseif ndims == 2
        plot_distance_2d(d; kwargs...)
    end
end

plot_distance_function(d::MobiusEuclidean; kwargs...) = plot_distance_2d(d; kwargs...)

# ------------------------------- connectivity ------------------------------- #
""" Plot connectivity matrix. """
function show_connectivity end

"""
    show_connectivity(W::Matrix)

Plot an entire connectivity matrix as a heatmap
"""
show_connectivity(W::Matrix; label=nothing, xlabel="neuron #", 
ylabel="neuron #", kwargs...) = heatmap(
    W, 
    xlabel=xlabel, 
    ylabel=ylabel, 
    aspect_ratio=:equal,
    colorbar=nothing,    
    label=label;
    kwargs...
)

show_connectivity!(W::Matrix; label=nothing) = heatmap!(
    W, 
    xlabel="neuron #", 
    ylabel="neuron #", 
    aspect_ratio=:equal,
    colorbar=nothing,    
    label=label,
    alpha=.33
)

""" 
    show_connectivity(W::Vector) 
Plot 1D connectivity matrix as a vector
"""
show_connectivity(W::Vector; label=nothing) = plot(
    W, 
    lw = 3,
    xlabel="neuron #", 
    ylabel="weights", 
    label=label
)

show_connectivity!(W::Vector; label=nothing) = plot!(
    W, 
    lw = 3, 
    xlabel="neuron #", 
    ylabel="weights", 
    label=label
)



""" 
    show_connectivity(W::Matrix, n::NTuple{N,Int}, i::Int) where N

show connectivity for a single neuron in a 2D lattice 
"""
function show_connectivity(W::Matrix, n::NTuple{N,Int}, i::Int) where N
    weights = reshape(W[i, :], n...)
    show_connectivity(weights)
end

"""
    show_connectivity(can::CAN, i::Int)

Show connectivity for a can's neuron given its index
"""
function show_connectivity(can::CAN, i::Int)
    if length(can.n) == 1
        p = plot()
        for (n, W) in enumerate(can.Ws)
            weights = reshape(W[i, :], can.n...)
            show_connectivity!(weights; label=nothing)
        end
        x = can.I[i]
        vline!([x[1]], label="Neuron $i", lw=4, color=:black)
    else
        p = plot()
        offsets = [[0, 0], [1, 0], [0, 1], [1, 1]]
        for (n, W) in enumerate(can.Ws)
            # plot connectivity map
            w = reshape(W[:, i], can.n...)

            Δx, Δy = can.n[2] * offsets[n][1], can.n[1] * offsets[n][2]
            x = collect(1:can.n[2]) .+ Δx
            y = collect(1:can.n[1]) .+ Δy
            heatmap!(x, y, w, 
                colorbar=nothing, 
                xaxis=false, 
                yaxis=false, 
                aspect_ratio=:equal,
                xticks=[], yticks=[]
            )
        end
        # separate heatmaps
        vline!([can.n[2]], lw=4, color=:white, label=nothing)
        hline!([can.n[1]], lw=4, color=:white, label=nothing)

        # mark the neruon's location
        for n in 1:length(can.Ws)    
            x = can.I[i]
            # Δ = reverse(can.n .* offsets[n])
            Δ = can.n .* offsets[n]
            scatter!(
                reverse(map(z->[z], x .+ Δ))..., 
                color=:green,
                label= nothing,
                ms=8)
        end

    end
    p
end

"""
    show_connectivity(can::CAN)

Show the connectivity for a can, the kernel and a 
few randomly selected neurons.
"""
function show_connectivity(can::CAN; kwargs...)
    idxs = rand(1:*(can.n...), 5)
    p = plot(can.kernel; title="Connectivity kernel")
    ps = map(i -> show_connectivity(can, i), idxs)
    plot(p, ps...;  size=(800, 600), kwargs...)
end
