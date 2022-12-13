import MyterialColors: red, green, black, white

# ---------------------------------- kernel ---------------------------------- #
function Plots.plot(K::AbstractKernel; σ = 4, kwargs...)
    x = -σ:0.001:σ |> collect
    y = K.(x)
    x̂ = x[argmin(y)]
    Plots.plot(
        x,
        y,
        xlabel = "distance Δx",
        ylabel = "w",
        lw = 2,
        label = nothing,
        grid = false,
        color = "black";
        kwargs...,
        # xlim = [-σ*abs(x̂), σ*abs(x̂)],
    )
end



# ---------------------------- distance functions ---------------------------- #
function plot_distance_1d(d::PeriodicEuclidean; kwargs...)
    x = 0:0.1:d.periods[1] |> collect

    p = plot(; xlabel = "x", ylabel = "distance", kwargs...)
    colors = [:black, :red, :green, :blue]
    for (i, point) in enumerate(rand(x, 4))
        y = [evaluate(d, point, xx) for xx in x]
        plot!(x, y, lw = 3, color = colors[i], label = nothing)
        vline!([point], color = colors[i], label = nothing)
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
    points = nothing,
    xlabel = "x",
    ylabel = "y",
    kwargs...,
)

    # create lattice of points
    X = (x × y) |> collect
    X = [[x...] for x in vec(X)]

    # get distance from a point
    points = isnothing(points) ? [[0, 0], rand(X, 3)...] : points
    pts = []
    for p in points
        Δx = [evaluate(d, p, x) for x in X]

        _plot = contourf(
            x,
            y,
            reshape(Δx, length(x), length(y))',
            aspect_ratio = :equal,
            xlabel = xlabel,
            ylabel = ylabel,
            linewidth = 0.25,
            grid = false,
        )
        scatter!([p[1]], [p[2]], color = :green, label = nothing, ms = 10, alpha = 1)
        push!(pts, _plot)
    end

    plt = plot(pts..., size = (600, 600); kwargs...)
    display(plt)
end

"""
    plot_distance_2d(d::PeriodicEuclidean; kwargs...)

Plot metrics of type PeriodicEuclidean
"""
function plot_distance_2d(d::PeriodicEuclidean; kwargs...)
    upperbound(x) = isfinite(x) ? x : 1
    # get coordinates mesh
    x = range(0, upperbound(d.periods[1]), length = 100) |> collect
    y = range(0, upperbound(d.periods[2]), length = 100) |> collect

    plot_distance_2d(d, x, y; kwargs...)
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


"""
plot_distance_function(d::MobiusEuclidean; kwargs...)

Plot distance for metric of type MobiusEuclidean
"""
function plot_distance_function(d::MobiusEuclidean; kwargs...)
    x = -1/2:0.075:1.2 |> collect
    y = 0:0.075:2π |> collect
    X = (x × y) |> collect
    X = [[x...] for x in vec(X)]

    plot_distance_2d(
        d,
        x,
        y;
        points = [[-1 / 2, 0], [0, 3], [0.5, 0.2], [0, 2π]],
        kwargs...,
    )
end

"""
plot_distance_function(d::SphericalAngle; kwargs...)

Plot distance for metric of type SphericalAngle
https://github.com/JuliaStats/Distances.jl/blob/master/src/haversine.jl

For an embedding of the sphere see: https://stackoverflow.com/questions/10473852/convert-latitude-and-longitude-to-point-in-3d-space
"""
function plot_distance_function(d::Union{SphericalDistance,SphericalAngle}; kwargs...)
    long = range(-π + 0.01, π - 0.01, length = 100) |> collect
    lat = range(-π / 2 + 0.01, π / 2 - 0.01, length = 100) |> collect
    X = (long × lat) |> collect
    X = [[x...] for x in vec(X)]

    points = [[-π, 0], [2, π / 2], [π - 1, π / 2 - 1], [π, -1]]
    plot_distance_2d(
        d,
        long,
        lat;
        points = points,
        xlabel = "longitude",
        ylabel = "latitude",
        kwargs...,
    )
end

# ------------------------------- connectivity ------------------------------- #
""" Plot connectivity matrix. """
function show_connectivity end

"""
    show_connectivity(W::Matrix)

Plot an entire connectivity matrix as a heatmap
"""
show_connectivity(
    W::Matrix;
    label = nothing,
    xlabel = "neuron #",
    ylabel = "neuron #",
    kwargs...,
) = heatmap(
    W,
    xlabel = xlabel,
    ylabel = ylabel,
    aspect_ratio = :equal,
    colorbar = nothing,
    label = label;
    kwargs...,
)

show_connectivity!(W::Matrix; label = nothing, kwargs...) = heatmap!(
    W,
    xlabel = "neuron #",
    ylabel = "neuron #",
    aspect_ratio = :equal,
    colorbar = nothing,
    label = label,
    alpha = 0.33;
    kwargs...,
)

""" 
    show_connectivity(W::Vector) 
Plot 1D connectivity matrix as a vector
"""
show_connectivity(W::Vector; label = nothing) =
    plot(W, lw = 3, xlabel = "neuron #", ylabel = "weights", label = label)

show_connectivity!(W::Vector; label = nothing) =
    plot!(W, lw = 3, xlabel = "neuron #", ylabel = "weights", label = label)



""" 
    show_connectivity(W::Matrix, n::NTuple{N,Int}, i::Int) where N

show connectivity for a single neuron in a 2D lattice 
"""
function show_connectivity(W::Matrix, n::NTuple{N,Int}, i::Int; kwargs...) where {N}
    weights = reshape(W[i, :], n...)
    show_connectivity(weights; kwargs...)
end

"""
    show_connectivity(can::CAN, i::Int)

Show connectivity for a can's neuron given its index
"""
function show_connectivity(can::CAN, i::Int; aspect_ratio = :equal, kwargs...)
    if can.d == 1
        p = plot(; kwargs...)
        for (n, W) in enumerate(can.Ws)
            weights = reshape(W[i, :], can.n...)
            show_connectivity!(weights; label = nothing, kwargs...)
        end
        x = can.I[i]
        vline!([x[1]], label = "Neuron $i", lw = 4, color = :black)
    elseif can.d == 2
        p = plot()
        offsets = map(offset_for_visual, can.offsets)
        for (n, W) in enumerate(can.Ws)
            # plot connectivity map
            w = reshape(W[:, i], can.n...)'

            Δx, Δy = can.n[2] * offsets[n][1], can.n[1] * offsets[n][2]
            x = collect(1:can.n[1]) .+ Δx
            y = collect(1:can.n[2]) .+ Δy
            heatmap!(
                x,
                y,
                w,
                colorbar = nothing,
                xaxis = false,
                yaxis = false,
                aspect_ratio = aspect_ratio,
                xticks = [],
                yticks = [];
                kwargs...,
                color = :bwr,
            )

        end

        # mark the neuron's location
        for n = 1:length(can.Ws)
            x = can.I[i]
            # Δ = reverse(can.n .* offsets[n])
            Δ = can.n .* offsets[n]
            scatter!(
                # reverse(map(z->[z], x .+ Δ))..., 
                map(z -> [z], x .+ Δ)...,
                color = :black,
                label = nothing,
                ms = 4,
            )
        end
    else
        error("Not implemented for d>2")
    end
    p
end

"""
    show_connectivity(can::CAN)

Show the connectivity for a can, the kernel and a 
few randomly selected neurons.
"""
function show_connectivity(
    can::CAN;
    idxs = nothing,
    aspect_ratio = :equal,
    size = (800, 600),
    kwargs...,
)
    idxs = isnothing(idxs) ? [1, rand(1:*(can.n...), 4)...] : idxs
    p = plot(can.kernel; title = "Connectivity kernel")

    clims = (minimum(can.Ws[1]), max(abs(minimum(can.Ws[1])) / 2, maximum(can.Ws[1])))
    ps = map(
        i -> show_connectivity(can, i; clims = clims, aspect_ratio = aspect_ratio),
        idxs,
    )
    plot(p, ps...; size = size, kwargs...) |> display
end



# ---------------------------------------------------------------------------- #
#                                   ONE FORMS                                  #
# ---------------------------------------------------------------------------- #
# -------------------------------- one oneform ------------------------------- #
function show_oneforms!(plt, ω::OneForm, xmin::Vector, xmax::Vector; dx = 3, x₀ = 0, y₀ = 0)
    length(xmin) != 2 && error("not implemented for d != 2")

    for x = xmin[1]:dx:xmax[1], y = xmin[2]:dx:xmax[2]
        o = ω([x, y])
        scatter!([x + x₀], [y + y₀], ms = 3.5, color = :black, label = nothing)
        plot!(
            plt,
            [x + x₀, x + x₀ + o[1]],
            [y + y₀, y + y₀ + o[2]],
            lw = 3,
            color = :black,
            label = nothing,
        )
    end

end


function show_oneforms(ω::OneForm, xmin::Vector, xmax::Vector; dx = 3)
    plt = plot(aspect_ratio = :equal, grid = false)
    show_oneforms!(plt, ω, xmin, xmax; dx = dx)
    plt
end


# ---------------------------- all ωᵢ for one CAN ---------------------------- #
function show_oneforms(can::CAN; kwargs...)
    plt = plot(; aspect_ratio = :equal, grid = false, size = can.n .* 10)

    x̄ = range(0, maximum(can.X[1, :]), length = can.n[1])
    ȳ = range(0, maximum(can.X[2, :]), length = can.n[2])
    xmin = [minimum(x̄), minimum(ȳ)]
    xmax = [maximum(x̄), maximum(ȳ)]

    for i = 1:can.d*2
        # get offset position for plotting
        offset = can.offsets[i] .* vec(maximum(can.X; dims = 2))

        # show one forms
        show_oneforms!(plt, can.Ω[i], xmin, xmax; x₀ = offset[1], y₀ = offset[2])
    end
    plt
end


# ----------------------- ωᵢ  over M (pullback from N) ----------------------- #
function show_oneforms(ω::OneForm, C::CoverSpace, args...; kwargs...)
    plt = plot(; aspect_ratio = :equal, grid = false, size = (800, 800))

    show_oneforms!(plt, ω, C, args...; kwargs...)
end


function show_oneforms!(
    plt,
    ω::OneForm,
    C::CoverSpace,
    xmin::Vector,
    xmax::Vector;
    dx = 10,
    scale = 1,
    color = :black,
    kwargs...,
)
    for x = xmin[1]:dx:xmax[1], y = xmin[2]:dx:xmax[2]
        # get ω at the corresponding mfld
        x̂, ŷ = C.ρ(x, y)
        o = ω([x̂, ŷ]) .* scale

        # vis
        scatter!(
            plt,
            [x],
            [y],
            ms = 3,
            msa = 0,
            msw = 0,
            markercolor = color,
            label = nothing;
            kwargs...,
        )
        plot!(
            plt,
            [x, x + o[1]],
            [y, y + o[2]],
            lw = 3,
            color = color,
            label = nothing;
            kwargs...,
        )
    end
    return plt
end



function plot_can_vector_fields!(plt, can, vel, x_actual, x_decoded)
    if can.d ==1
        x0, x1 = minimum(can.X), maximum(can.X)
        for x in range(x0, x1, length=50)
            y = can.Ω[1]([x])[1][1] * 4
            plot!(plt, [x, x], [0, y], 
                ylim=[0, 2],
                lw=2, color=:red, alpha=.5, label=nothing)
        end

        plot!(plt, [x_actual[1], x_actual[1]], [0, 1.5],
            lw=4, color=:black, label=nothing
        )
        plot!(plt, [x_decoded[1], x_decoded[1]], [0, 1.5],
            lw=4, color=:red, label=nothing
        )
        return
    end
    n = size(can.X, 2)
    scaling = 25.0
    
    colors = [white, white, red, red, green, green]
    for i in 1:4:n
        x = can.X[:, i]
        scatter!(plt, [[x] for x in x]..., label=nothing, color=:black, ms=3)

        for (j, o) in enumerate(can.Ω)
            j % 2 == 0 && continue
            v = o(x) * vel[(Int ∘ ceil)(j/2)]
            plot!(plt,
                [x[1], x[1]+v[1]*scaling],
                [x[2], x[2]+v[2]*scaling],
                lw=2, color=colors[j], label=nothing
            )
        end
    end
end