using Plots


# ---------------------------------- kernel ---------------------------------- #
function Plots.plot(K::Kernel)
    x = -100:.1:100 |> collect
    y = K.(x)
    Plots.plot(x, y, xlabel="distance Δx", ylabel="w", lw=2, color="black")
end


# ------------------------------- connectivity ------------------------------- #
""" Plot connectivity matrix. """
function show_connectivity end

"""
    show_connectivity(W::Matrix)

Plot an entire connectivity matrix as a heatmap
"""
show_connectivity(W::Matrix) = heatmap(W, xlabel="neuron #", ylabel="neuron #", aspect_ratio=:equal, size=(500, 500))

function show_connectivity(W::Matrix, n::NTuple{N,Int}, i::Int) where N
    weights = reshape(W[i, :], n...)
    p = show_connectivity(weights)
end


function show_connectivity(can::CAN, i::Int)
    weights = reshape(can.W[i, :], can.n...)
    p = show_connectivity(weights)
    x = can.I[i]
    scatter!(reverse(map(z->[z], x))..., color=:green, label="p", ms=12)
    p
end

function show_connectivity(can::CAN)
    idxs = rand(1:*(can.n...), 6)
    ps = map(i -> show_connectivity(can, i), idxs)
    plot(ps...)
end
