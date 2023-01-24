using Term.Progress
include("settings.jl")

import MyterialColors: blue_grey_dark, salmon
using NearestNeighbors

"""
Try to estimate the dimensionality of the data and the intrinsic
dimensionality of the activity manifold with PCA/local PCA. 

PLOT_ACTIVITY_PCA_VARIANCE_EXPLAINED
Plot the fraction of variance explained by the first 250PCs for each 
network's activity.


ESTIMATE_LOCAL_PCA_PARAMS_SENSITIVITY
Estimate the sensitivity of the local PCA approach to intrinsic dimensionality
to parameters by doing a paramter sweep on torus simulation data.


ESTIMATE_INTRINSIC_DIMENSIONALITY
Run analysis on (10d embedded) data for each network.
"""

import GeneralAttractors.Analysis.ManifoldAnalysis: 
        fraction_variance_explained, find_fraction_variance_explained_elbow, pca_dimensionality_reduction
move_to_datadir(supervisor, "mfld_top")

PLOT_ACTIVITY_PCA_VARIANCE_EXPLAINED = false
ESTIMATE_LOCAL_PCA_PARAMS_SENSITIVITY = false
ESTIMATE_INTRINSIC_DIMENSIONALITY = true

if PLOT_ACTIVITY_PCA_VARIANCE_EXPLAINED
    dim_est_params = AnalysisParameters(
        max_nPC = 250,  # max num of PCs
        pca_pratio = 0.999999,       # fraction of variance explained
    )


    plt = plot(
        xlabel = "PC", ylabel = "Fraction of variance explained",
        grid = false,
        size = (1000, 600)
    )

    for (network, color) in zip(networks, networks_colors)
        print(hLine(network; style="red"))
        filters = Dict{Symbol, Any}(
            :tag => "manifold_topology_data",
            :can => network,
        )

        X = load_and_concat_activations(; filters...) 
        pca = pca_dimensionality_reduction(X, dim_est_params)[1]

        # get fraction of variance explained and dimensionality
        fvariance_explained = fraction_variance_explained(pca) 
        # plot
        σ = cumsum(fvariance_explained)

        plot!(
            1:length(σ), σ,
            lw = 1.5,
            label = nothing,
            color = color,
        )
        scatter!(
            1:length(σ), σ,
            markersize = 3.0,
            label = network,
            color = color,
            msa=0, msw=0,
        )
        print(hLine(; style="dim"))
    end

    save_plot(supervisor, plt, "04c_mfld_top_estimate_dimensionality")
end


if ESTIMATE_LOCAL_PCA_PARAMS_SENSITIVITY
    network = "torus"
    # load low dimensional embedding data
    filters = Dict{Symbol, Any}(
        :tag => "d10_embeddings",
        :can => network,
    )
    _, M = ProjectSupervisor.fetch(supervisor; filters...) 
    @assert length(M) == 1 length(M)
    M = M[1]
    @assert size(M, 1) == 10 size(M)
    @info "Loaded $(size(M, 2)) samples of dim. $(size(M, 1))"

    # build NN tree
    nntree = KDTree(M; reorder = false, leafsize = 5)

    # params
    neighborhood_sizes = 25:25:1000
    pratios = (.5, .75, .9,  .99)
    colors = getfield.(Palette(blue_grey_dark, salmon; N = length(pratios)).colors, :string)

    plt = plot(
        xlabel = "Neighborhood size", ylabel = "Fraction of variance explained",
        grid = false,
        size = (1000, 600),
        title = network,
    )
    for (p, color) in zip(pratios, colors)
        estimated_dim = []
        for ns in neighborhood_sizes
            prms = AnalysisParameters(
                intrinsic_d_nseeds = 250,
                intrinsic_d_pratio = p,
                intrinsic_d_neighborhood_size = ns;
            )
            
            push!(
                estimated_dim, 
                estimate_intrinsic_dimensionality(M, prms; nntree=nntree) |> mean
            )
        end

        plot!(
            neighborhood_sizes,
            estimated_dim,
            lw = 3,
            label = "frac. var. expl. = $p",
            color = color,
        )

    end

    hline!(plt, [2], label = "True dimensionality", color = "black", lw = 3, alpha=.5, ls = :dash)
    save_plot(supervisor, plt, "04c_mfld_top_estimate_dimensionality_params_sensitivity")
    display(plt)
end



if ESTIMATE_INTRINSIC_DIMENSIONALITY
    for network in networks
        filters = Dict{Symbol, Any}(
            :tag => "d10_embeddings",
            :can => network,
        )
        _, M = ProjectSupervisor.fetch(supervisor; filters...) 
        @assert length(M) == 1 length(M)
        M = M[1]
        @assert size(M, 1) == 10 size(M)


        prms = AnalysisParameters(
                intrinsic_d_nseeds = 250,
                intrinsic_d_pratio = 0.75,
                intrinsic_d_neighborhood_size = 500;
            )
        d = estimate_intrinsic_dimensionality(M, prms)
        @info "Estimated intrinsic dimensionality of $network: $(d |> mean) ± $(d |> std)"
    end
end