using ManifoldLearning
using NearestNeighbors
using NPZ
using Statistics: mean
using MultivariateStats
using Plots

""" 
Preprocess neural data

Start by temporarily downsampling the data
through averaging, followed by PCA projection
on the first X PCs.
Return the projected data.
"""
function preprocessing()::Matrix
    S = npzread("S.npz");  # n neurons × n time steps
    @info "Loaded S" size(S)

    # take the average every 25 ms
    idxs = 1:5:size(S, 2) |> collect
    S̄ = hcat(map(i->mean(S[:, i:i+4], dims=2), idxs)...)
    S = 0  # remove memory usage

    # get first 10 PCs
    @info "Fitting PCA on reduced S" size(S̄)
    pca = fit(PCA, S̄, maxoutdim=5);

    # project with PCA
    S_pca = predict(pca, S̄ ) 
    npzwrite("S_pca_space.npz", S_pca)
    S_pca
end


""" further reduce dimensionality with isomap """
function embed3d()
    S = npzread("S_pca_space.npz")
    @info "embedding S" size(S)

    iso = fit(Isomap, S,  k=10, maxoutdim=3)
    M = predict(iso, S)
    @info "done M" size(M)
    npzwrite("M.npz", M)
    return M
end


# -------------------------------- processing -------------------------------- #
preprocessing()
embed3d()


# --------------------------------- visualize -------------------------------- #
# Vidx = npzread("Vidx.npz")[1:5:end][1:25:end]
# V = npzread("V.npz")[:, 1:5:end][:, 1:25:end]
# M̄ = npzread("M.npz")

# # plt = scatter3d(M̄[1, :],
# #         M̄[2, :],
# #         M̄[3, :],
# #         color=:black, label=nothing, alpha=.02,
# #         camera=(135, 15))

# # # plot traces for each velocity vector
# # for v in unique(V)
# #     idx = findall(V .== v)
# #     scatter3d!(M̄[1, idx],
# #             M̄[2, idx],
# #             M̄[3, idx],
# #             label=string(v)
# #             )
# # end

# # display(plt)  