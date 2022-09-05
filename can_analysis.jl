using ManifoldLearning
using NearestNeighbors
using NPZ
using Statistics: mean
using MultivariateStats

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

    # keep only neurons with a large enough mean activity
    μ = mean(S, dims=2)
    S = S[[(μ .> .005)...], :]

    # take the average every 25 ms
    idxs = 1:5:size(S, 2) |> collect
    S̄ = hcat(map(i->mean(S[:, i:i+4], dims=2), idxs)...)
    S = 0  # remove memory usage

    # get first 10 PCs
    @info "Fitting PCA on reduced S" size(S̄)
    pca = fit(PCA, S̄, maxoutdim=5);

    # project with PCA
    predict(pca, S̄ ) 
end


# S = preprocessing();
# npzwrite("S_pca_space.npz", S)

# further reduce dimensionality with isomap
S = npzread("S_pca_space.npz")
iso = fit(Isomap, S[:, 1:100:end],  k=10, maxoutdim=3)
M = predict(iso, S[:, 1:100:end])
npzwrite("M.npz", M)

# visualize
scatter3d(M[1, 1:100:end], M[2, 1:100:end], M[3, 1:100:end])