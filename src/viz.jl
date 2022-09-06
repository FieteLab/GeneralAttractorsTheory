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
show_connectivity(W::Matrix) = heatmap(W, xlabel="neuron #", ylabel="neuron #", aspect_ratio=:equal, size=(500, 500)) |> display



