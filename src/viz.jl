using Plots
using Measures

# ---------------------------------- kernel ---------------------------------- #
function Plots.plot(K::Kernel; kwargs...)
    x = -10:.1:10 |> collect
    y = K.(x)
    Plots.plot(x, y, xlabel="distance Δx", ylabel="w", lw=2, color="black"; kwargs...)
end


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
            show_connectivity!(weights; label="Offset $n")
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
                label= n == 1 ? "p: $x" : nothing,
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
function show_connectivity(can::CAN)
    idxs = rand(1:*(can.n...), 5)
    p = plot(can.kernel; title="Connectivity kernel")
    ps = map(i -> show_connectivity(can, i), idxs)
    plot(p, ps...;  size=(800, 600))
end
