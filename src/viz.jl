using Plots


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
show_connectivity(W::Matrix) = heatmap(
    W, 
    xlabel="neuron #", 
    ylabel="neuron #", 
    aspect_ratio=:equal,
    colorbar=nothing,    
)

""" 
    show_connectivity(W::Vector) 
Plot 1D connectivity matrix as a vector
"""
show_connectivity(W::Vector) = plot(
    W, 
    xlabel="neuron #", 
    ylabel="weights", 
    label=nothing
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
    weights = reshape(can.W[i, :], can.n...)
    p = show_connectivity(weights)

    if length(can.n) > 1
        x = can.I[i]
        scatter!(reverse(map(z->[z], x))..., color=:green, label="p: $x", ms=12)
    else
        x = can.I[i]
        vline!([x[1]], label="Neuron $i")
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
    plot(p, ps...;  size=(800, 800))
end
