include("settings.jl")


import MyterialColors: blue_grey_dark, salmon
using NearestNeighbors

"""
Try to estimate the dimensionality of the data and the intrinsic
dimensionality of the activity manifold with PCA/local PCA. 

PLOT_EXTRINSIC_DIMENSIONALITY
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
# move_to_datadir(supervisor, "mfld_top4")

PLOT_EXTRINSIC_DIMENSIONALITY = true
ESTIMATE_LOCAL_PCA_PARAMS_SENSITIVITY = false
ESTIMATE_INTRINSIC_DIMENSIONALITY = true



if PLOT_EXTRINSIC_DIMENSIONALITY
    dim_est_params = AnalysisParameters(
        max_nPC = nothing,  # Use all PCs
        pca_pratio = 1.0,  # Use all components
    )

    # Create individual plots for each network
    for (network, color) in zip(networks, networks_colors)
        print(hLine(network; style="red"))
        filters = Dict{Symbol, Any}(
            :tag => "manifold_topology_data",
            :can => network,
            :η => 0.0,
        )

        X = nothing
        try
            X = load_and_concat_activations(; filters...) 
        catch e
            @warn "Could not load data for network $network" e
            continue
        end
        pca = pca_dimensionality_reduction(X, dim_est_params)[1]
        
        # get fraction of variance explained for first 100 PCs
        fvariance_explained = fraction_variance_explained(pca)[1:25]

        plt = plot(
            xlabel = "PC", 
            ylabel = "Variance explained (%)",
            grid = false,
            size = (800, 500),
            title = network,
            ylims = (0, 100);
            xlims = (0, 25),
            plot_font_size_kwargs...
        )

        # Plot as bar plot
        plot!(
            plt,
            1:length(fvariance_explained),
            fvariance_explained,
            seriestype=:bar,
            lw = 0,
            label = nothing,
            color = color,
            alpha = 0.6
        )

        display(plt)
        save_plot(supervisor, plt, "f3_D_extrinsic_dimensionality_$(network)")
    end
end


if ESTIMATE_LOCAL_PCA_PARAMS_SENSITIVITY
    network = "torus"
    # load low dimensional embedding data
    filters = Dict{Symbol, Any}(
        :tag => "d10_embeddings",
        :can => network,
        :η => 0.0,
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
        title = network;
        plot_font_size_kwargs...
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
    save_plot(supervisor, plt, "f3_D_intrisic_dim_params_sensitivity")
    display(plt)
end



if ESTIMATE_INTRINSIC_DIMENSIONALITY
    # Create individual plots for each network
    for (network, color) in zip(networks, networks_colors)
        println("Plotting $network")
        filters = Dict{Symbol, Any}(
            :tag => "d10_embeddings",
            :can => network,
            :η => 0.0,
            :extension => "npz",
        )

        df, M = ProjectSupervisor.fetch(supervisor; filters...) 
        @assert length(M) == 1 length(M)
        M = M[1]
        @assert size(M, 1) == 10 size(M)

        # Get dimensionality and variances
        d, all_variances = estimate_intrinsic_dimensionality(M, intrinsic_dimensionality_prms)
        @info "Estimated intrinsic dimensionality of $network: $(d |> mean) ± $(d |> std)"

        # Calculate mean and std of variance explained for all PCs
        n_dims = size(M, 1)  # Use full dimensionality
        var_means = zeros(n_dims)
        var_stds = zeros(n_dims)
        
        for pc in 1:n_dims
            pc_vars = [v[pc] for v in all_variances if length(v) >= pc]
            var_means[pc] = mean(pc_vars)
            var_stds[pc] = std(pc_vars)
        end

        # Create individual plot
        plt = plot(
            xlabel = "PC",
            ylabel = "Variance explained (%)",
            grid = false,
            size = (800, 500),
            title = network,
            ylims = (0, 100),
            xlims = (0, 25),
        )

        # Plot all PCs
        plot!(
            plt,
            1:n_dims,
            var_means,
            yerror=var_stds,
            color=color,
            label=nothing,
            seriestype=:bar,
            alpha=0.6
        )

        display(plt)
        save_plot(supervisor, plt, "f3_D_intrisic_dim_local_pca_$(network)")
    end
end