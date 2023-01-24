using GeneralAttractors
using Distances
import GeneralAttractors: plot_distance_function, SphericalDistance, MobiusEuclidean

include("settings.jl")
using LinearAlgebra

"""
Plot eigenspectrum of connectivity matrices.
"""

# ----------------------------- plot connectivity ---------------------------- #
plots = []

for network in networks 
    can = network_makers[network](:single; n=64)

    w = reshape(can.W[1,:], can.n)' |> Matrix
    # w = can.W

    r = rank(w)
    println(network, " - rank: ", r)

    evals = eigvals(w)
    # evecs = eigvecs(w)


    x = real.(evals)
    y = imag.(evals)

    plt = scatter(
        x, y,
        color=:red,
        ms=2,
        msw=0,
        msa=0,
        label=nothing,
        title= network * " - rank: $r",
        aspect_ratio=:equal,
    )

    plot!(
        sin.(range(0, 2π; length=100)),
        cos.(range(0, 2π; length=100)),
        color=:black,
        linewidth=1,
        label=nothing,
        alpha=.5
    )

    push!(plots,    
        plt,
    )
 
    end

fig = plot(plots..., size=(500, 500))
fig 