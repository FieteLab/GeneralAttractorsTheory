
"""
Analysis of manifold topology and dimensionality
from point cloud data (e.g. neural activity recordings)
"""
module ManifoldAnalysis
using ManifoldLearning
using NearestNeighbors
import NearestNeighbors: NNTree
using MultivariateStats
using Plots
using Term.Progress
using Term
using Ripserer
using PersistenceDiagrams: PersistenceDiagram
using Statistics


import ...Simulations: History
import ...Analysis: AnalysisParameters

include("analysis_viz.jl")

export pca_dimensionality_reduction,
    isomap_dimensionality_reduction,
    estimate_intrinsic_dimensionality

# ----------------------------------- utils ---------------------------------- #
""" compute fraction of variance explained by each PC """
fraction_variance_explained(M::PCA) = principalvars(M) ./ tvar(M) * 100

"""
    find_fraction_variance_explained_elbow(σ::Vector{Float64})::Int

Given a vector σ with the fraction of variance explained of each PC, 
find the "elbow" to identify the "relevant" number of PC.

See: https://ieeexplore.ieee.org/document/5961514
https://towardsdatascience.com/detecting-knee-elbow-points-in-a-graph-d13fc517a63c
"""
function find_fraction_variance_explained_elbow(σ::Vector{Float64})::Int
    length(σ) == 1 && return 1
    # check input
    @assert σ[1] > σ[end] "Fraction of variance explained should be decreasing, did you `cumsum`?"

    # fit a line through the first and last point
    # consider that σ₁ = f(x₁) - x₁=1
    # we don't start at x=0 so fix the intercept
    m = (σ[end] - σ[1]) / (length(σ) - 1)
    y₀ = σ[1]
    ŷ(x) = m * x + y₀

    # compute the line between the start and end value
    x = 0:length(σ)-1
    y = ŷ.(x)

    # get the vertical distance between σ and the line
    Δ = σ .- y

    # get the point with the largest vertical distance
    return argmax(Δ)
end


"""
Take the average/sum of the population activity across
copies of the network
"""
function population_average(history::History; skip = 10)::Matrix
    n, _, m = size(history.S)
    S = real.(reshape(mean(history.S, dims = 2), (n, m)))
    not_nan_cols = map(c -> !any(isnan.(c)), eachcol(S)) |> collect
    S = S[:, not_nan_cols]
    S = S[:, skip:end-skip]
    return S
end

# ---------------------------------------------------------------------------- #
#                                      PCA                                     #
# ---------------------------------------------------------------------------- #




function pca_dimensionality_reduction(
    S::Matrix,
    params::AnalysisParameters = AnalysisParameters();
)
    # get number of PCs
    nPC = if !isnothing(params.max_nPC)
        params.max_nPC
    else
        25
    end

    # reduce dimensionality with PCA
    @info "PCA dimensionality reduction" size(S) params.max_nPC params.pca_pratio nPC
    pca_model = fit(PCA, S; pratio = params.pca_pratio, maxoutdim = nPC)
    S_pca_space = predict(pca_model, S)
    @info "pca fitting completed $(length(principalvars(pca_model))) PCs" size(S_pca_space)

   return pca_model, S_pca_space
end


# ---------------------------------------------------------------------------- #
#                                    ISOMAP                                    #
# ---------------------------------------------------------------------------- #

"""
    function isomap_dimensionality_reduction(
        X::Matrix,
        params::AnalysisParameters = AnalysisParameters();
        visualize = true,
    )::Nothing

Embed a data matrix X of shape n_cells × n_timepoints
into a lower dimensional space using ISOMAP.
"""
function isomap_dimensionality_reduction(
    X::Matrix,
    params::AnalysisParameters = AnalysisParameters();)
    # fit
    @info "Performing ISOMAP" size(X) params.n_isomap_dimensions params.isomap_downsample
    iso = ManifoldLearning.fit(
        Isomap,
        X[:, 1:params.isomap_downsample:end];
        k = params.isomap_k,
        maxoutdim = params.n_isomap_dimensions,
    )
    M = predict(iso, X)
    @info "isomap fitting completed" size(M)


    return iso, M
end



# ---------------------------------------------------------------------------- #
#                                      TDA                                     #
# ---------------------------------------------------------------------------- #
"""
        function run_tda(
            simulation_name::String, 
            params::AnalysisParameters=AnalysisParameters(),
        )    

    Run Topological data analysis on simulation data in PCA space. 
    This is used to reconstruct the topology of the activity manifold. 
    It might be very slow for large number of dimensions.

    Ref: https://mtsch.github.io/Ripserer.jl/dev/
"""
function tda_on_pointcloud(X::Matrix, params::AnalysisParameters)
    # convert M in a vector of tuples for TDA
    n = (Int ∘ round)(size(X, 2) / params.tda_downsample_factor)
    X̄ = [Tuple(x) for x in rand(collect(eachcol(X)), n)]

    # fit TDA
    @info "Fitting TDA on X̄" length(X̄) length(X̄[1])
    tda = ripserer(
        X̄;
        dim_max = params.tda_dim_max,
        verbose = true,
        reps = false,
        threshold = params.tda_threshold,
    )

    # plot results
    plt = plot(plot(tda), barcode(tda), size = (1000, 800))
    return tda, plt
end




# ---------------------------------------------------------------------------- #
#                                   LOCAL PCA                                  #
# ---------------------------------------------------------------------------- #

"""
function estimate_intrinsic_dimensionality(
    simulation_name::String,
    params::AnalysisParameters=AnalysisParameters();
    verbose::Bool = true
)::Vector{Int}

Fit PCA to local neighborhoods on the data manifold to estimate
intrinsic dimensionality. 
"""
function estimate_intrinsic_dimensionality(
    M::Matrix,
    params::AnalysisParameters = AnalysisParameters();
)::Vector{Int}
    @info "Estimating intrinsic dimensionality" size(M) params.intrinsic_d_nseeds params.intrinsic_d_neighborhood_size

    # build nearest neighbor tree
    nntree = KDTree(M; reorder = false, leafsize = 5)

    # sample random points on the manifold
    seeds_idxs = rand(1:size(M, 2), params.intrinsic_d_nseeds)
    seeds = M[:, seeds_idxs]

    # get neighborhoods
    k = (Int ∘ round)(params.intrinsic_d_neighborhood_size)
    Us, _ = knn(nntree, seeds, k)

    # for each neighborhood fit PCA and get the number of PCs
    D = []  # store "dimensionality" of each tangent vector space
    for U in Us
        @assert length(U) == k
        pca_model = fit(PCA, M[:, U]; pratio = 0.999999, maxoutdim = size(M, 2))
        d = find_fraction_variance_explained_elbow(principalvars(pca_model))
        push!(D, d)
    end
    return D
end


end
