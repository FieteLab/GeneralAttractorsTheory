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


function plot_distance_2d(d::PeriodicEuclidean; kwargs...)
    @info "distance" d.periods
    upperbound(x, ) = isfinite(x) ? x : 2
    # get coordinates mesh
    x = range(0, upperbound(d.periods[1]), length=100) |> collect
    y = range(0, upperbound(d.periods[2]), length=100) |> collect
    X = (x × y) |> collect
    X = [[x...] for x in vec(X)]

    # get distance from a point
    pts = []
    for _ in 1:4
        p = [rand(x), rand(y)]
        Δx = [evaluate(d, p, x) for x in X]
        Δx = reshape(Δx, length(x), length(y))'
        _plot = heatmap(x, y, Δx, aspect_ratio=:equal, grid=false,
            xlabel="x", ylabel="x"
        )
        scatter!([p[1]], [p[2]], color=:green, label=nothing, ms=10, alpha=.5)
        push!(pts, _plot)
    end

    plt = plot(pts..., size=(600, 600); kwargs...)
    display(plt)
end


function plot_distance_2d(d::MobiusEuclidean; kwargs...)
    x = range(0, 2π, length=100) |> collect
    y = range(0, 1, length=100) |> collect
    X = (x × y) |> collect
    X = [[x...] for x in vec(X)]

    # get distance from a point
    pts = []
    points = [[0, 0], [.5, .2], [.05, .2], [2π-.05, .2]]
    for p in points
        # p = [rand(x), rand(y)]
        Δx = [evaluate(d, p, x) for x in X]
        Δx = reshape(Δx, length(x), length(y))'
        
        _plot = heatmap(
            x, y, Δx,  grid=false,
            xlabel="x", ylabel="x"
        )
        scatter!([p[1]], [p[2]], color=:green, label=nothing, ms=8, alpha=.8)
        push!(pts, _plot)
    end

    plt = plot(pts..., size=(600, 600); layout=(4,1), kwargs...)
    display(plt)
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
    p = show_connectivity(weights)
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
            w = reshape(W[i, :], can.n...)

            Δx, Δy = can.n[1] * offsets[n][1], can.n[2] * offsets[n][2]
            x = collect(1:can.n[1]) .+ Δx
            y = collect(1:can.n[2]) .+ Δy
            heatmap!(x, y, w, 
                colorbar=nothing, 
                xaxis=false, 
                yaxis=false, 
                aspect_ratio=:equal,
                xticks=[], yticks=[]
            )

            x = can.I[i]
            Δ = reverse(can.n .* offsets[n])
            scatter!(
                reverse(map(z->[z], x .+ Δ))..., 
                color=:green,
                label= n == 1 ? "neuron" : nothing,
                ms=12)
        end
        vline!([can.n[1]], lw=4, color=:white, label=nothing)
        hline!([can.n[1]], lw=4, color=:white, label=nothing)
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
