using Term.Progress
include("settings.jl")

"""
Fit a PCA model to a data from a large number of simulations to 
estimate dimensionality. 
"""

import GeneralAttractors.Analysis.ManifoldAnalysis: 
        fraction_variance_explained, find_fraction_variance_explained_elbow, pca_dimensionality_reduction
move_to_datadir(supervisor, "mfld_top")


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
plt
