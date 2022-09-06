"""
Continuous Attractor Network for grid cells. 

Based on Burak & Fiete 2009.

"""

import Base.Iterators: product as ×  # cartesian product
using Plots
using Term.Progress
import Term
import Distances: PeriodicEuclidean
import Distances
import LinearAlgebra: ⋅
using NPZ
using Statistics
import Base: Iterators


""" compute weight for an entry in the connectivity matrix """
W₀(x::Float64) = a * exp(-γ*abs(x)^2) - exp(-β*abs(x)^2)
W₀(x::Vector{Float64}) = W₀.(x) |> sum

function plot_connectivity(W, X)
    plots = []
    for i in (1, 15, 1500, 1921)
        weights = reshape(W[i, :], n, n)
        p = heatmap(weights,  colorbar=nothing, aspect_ratio=:equal, xlabel="neuron #", ylabel="neruon #", title= "neuron $i")

        x, y = X[i] .+ n/2 .+ 1
        scatter!([y], [x],  ms=10, color="green", label=nothing)
        push!(plots, p)
    end
    plot(plots..., size=(500, 500)) |> display
end


function get_Θ()
    # specify the orientation of each neuron (2x2 tiled square of orientations)
    orientations = [[[0.0, 1.0], [1.0, 0.0]] [[0.0, -1.0], [-1.0, 0.0]]] # 1=North, 2=East, 3=Sout, 4=West
    Θ = repeat(orientations, Int(n/2), Int(n/2))  # n × n mtx
    reshape(Θ, n^2)  # n^2 vector
end

function make_network()
    Θ = get_Θ()
    
    # get the position of every neuron in the lattice
    n̂ = Int(n/2)
    X::Matrix{Vector} = [[x...] for x in ((-n̂:(n̂-1)) × (-n̂:( n̂-1)))]

    # make recurrent weight matrix
    W::Matrix{Float64} = zeros(Float64, n^2, n^2)
    metric = PeriodicEuclidean([n, n])
    
    pbar = ProgressBar(columns=:detailed, expand=false)
    Term.Progress.with(pbar) do
        ijob = addjob!(pbar; description="Computing W", N=n^2)
        for i in 1:n^2
            for j in 1:n^2   
                # j > i && continue
                x = Distances.evaluate(metric, X[i] - (l .* Θ[i]), X[j])
                W[i, j] = W₀(x)
                # W[j, i] = W₀(x)
            end
            update!(ijob)
        end
    end
    plot_connectivity(W, X)
    npzwrite("W.npz", W)
    return Θ, W, X
end



function run()
    W = npzread("W.npz")
    Θ = get_Θ()
    S = zeros(Float64, n^2)  # initialize the state of all neurons
    T = 0:dt:150 |> collect
    
    pause = (false, 250, [0.0, 0.0])
    flow_directions = [
        (false, 150, [cos(0), sin(0)]),
        (false, 150, [cos(π/5), sin(π/5)]),
        (false, 150, [cos(π/2-π/5), sin(π/2-π/5)]),

        pause,
        (true, 10_000, [cos(0), 0.0] .* .6),
        pause,
        (true, 10_000, [0.0, sin(π/2)] .* .6),
        pause,
        (true, 10_000, [-cos(0), 0.0] .* .6),
        pause,
        (true, 10_000, [0.0, -sin(π/2)] .* .6),

        pause,
        (true, 10_000, [cos(π/4), sin(π/4)] .* .6),
        pause,
        (true, 10_000, [cos(-π/4), sin(-π/4)] .* .6),
        pause,
    ]

    # random = (true, 2000, nothing)
    # randos = vcat(repeat([random, pause], 250)...)
    # flow_directions = vcat(flow_directions..., randos...)

    T = sum(map(x->x[2], flow_directions)) / dt
    Ts = sum(
        map(x -> x[1]*x[2], flow_directions)
    ) / dt |> Int

    # initialize data storage
    storage = zeros(Float64, length(S), Ts)
    v_storage = zeros(Float64, 2, Ts)
    vidx_storage = zeros(Int, Ts)

    pbar = ProgressBar()
    anim = Animation()
    Term.Progress.with(pbar) do
        job = addjob!(pbar, description="Simulation",  N=Int(T))
        framen, saved_frames = 0, 0
        for (vv, (r, _t, v)) in enumerate(flow_directions)
            v = if isnothing(v)
                θ_random = rand(0:.5:2π)
                scale = rand(.2:.1:.8)
                [cos(θ_random), sin(θ_random)] .* scale
            else
                v
            end

            B = 1 .+ α.*map(θ -> θ ⋅ v, Θ) 
            for t in 1:dt:_t    
                # get random input
                η = rand(Float64, length(S)) * 0.05

                # compute activations
                S += (f.(W * S + η + B) - S) / τ
                update!(job)

                # store data
                r && begin
                    saved_frames += 1
                    storage[:, saved_frames] = S
                    v_storage[:, saved_frames] = v
                    vidx_storage[saved_frames] = vv
                end

                # animation
                # framen += 1
                # elapsed = dt * framen
                # if elapsed % 25 == 0
                #     ttl = (round(v[1]; digits=2), round(v[2]; digits=2))
                #     heatmap(
                #         reshape(S, n, n), 
                #         colorbar=nothing, title="Time: $elapsed ms | v $vv:$ttl",
                #         clims=(0.0, 0.5),
                #         aspect_ratio=:equal, size=(400, 400)
                #     )
                #     scatter!([5], [5], ms=25, color="white", alpha=.5, label=nothing)
                #     plot!([5, 5+(v[2]* n/6)], [5, 5+(v[1]* n/6)], lw=8, color="green", label="v⃗")
                #     scatter!([5], [5], ms=10, color="white", label=nothing)
                #     frame(anim)
                # end
            end
        end

    end

    npzwrite("S.npz", storage)
    npzwrite("V.npz", v_storage)
    npzwrite("Vidx.npz", vidx_storage)
    # gif(anim, "test2.gif", fps=20)
end



# --------------------------------- analysis --------------------------------- #
function plot_activity_flattened_over_time(S, V, Vidx)
    p1 = heatmap(Vidx[1:25:length(Vidx), :], colorbar=nothing, 
            ylabel="frame #", xticks=nothing,
            c=:grays, title="v⃗ id"
    )

    p2 = heatmap(V[:, 1:25:length(Vidx)]', colorbar=nothing, 
                yticks=nothing, xticks=nothing,
                clim=(-1, 1), c=:bwr, title="v⃗"
    )

    p3 = heatmap(S[1:12n, 1:25:length(Vidx)]', colorbar=nothing, 
                xlabel="cell #", yticks=nothing, title="activity"
    )

    plot(p1, p2, p3, size=(1200, 700),
        layout = grid(1, 3, widths=[.1, .1, .8])
    ) |> display

end


# ---------------------------------------------------------------------------- #
#                                      RUN                                     #
# ---------------------------------------------------------------------------- #

# parameters
n = 64  # √ of number of neurons
f(x) = max(0, x)  # relu
τ = 10  # ms
dt = 0.5 # ms

a = 1.0
α = 0.10315
l = 1.0
λ = 13
β = 3/(λ^2)
γ = 1.05 * β



# create network and simulate
Θ, W, X = make_network()
run()

