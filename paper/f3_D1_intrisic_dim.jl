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
move_to_datadir(supervisor, "mfld_top4")

PLOT_EXTRINSIC_DIMENSIONALITY = false
ESTIMATE_LOCAL_PCA_PARAMS_SENSITIVITY = false
ESTIMATE_INTRINSIC_DIMENSIONALITY = true



if PLOT_EXTRINSIC_DIMENSIONALITY
    dim_est_params = AnalysisParameters(
        max_nPC = 150,  # max num of PCs
        pca_pratio = 0.999999,       # fraction of variance explained
    )


    plt = plot(
        xlabel = "PC", ylabel = "Fraction of variance explained",
        grid = false,
        size = (1000, 600);
        plot_font_size_kwargs...
    )

    for (network, color) in zip(networks, networks_colors)
        # network != "cylinder" && continue
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
        # scatter!(
        #     1:length(σ), σ,
        #     markersize = 3.0,
        #     label = network,
        #     color = color,
        #     msa=0, msw=0,
        # )

        above = findfirst(σ .> 80)
        println(above, network)

        scatter!(
            [above], [σ[above]],
            markersize = 5.0,
            label = nothing,
            color = color,
            msa=0, msw=0,
        )
        plot!(
            [above, above], [0, σ[above]],
            lw = 1.5,
            label = nothing,
            color = color,
            ls = :dash,
        )
        annotate!((above, 0, network))


        print(hLine(; style="dim"))
    end

    display(plt)
    save_plot(supervisor, plt, "f3_D_extrinsic_dimensionality")
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
    # Create a figure with subplots for each network
    n_networks = length(networks) - 1  # excluding line
    n_cols = 2
    n_rows = ceil(Int, n_networks/n_cols)
    plt = plot(
        layout=(n_rows, n_cols),
        size=(1000, 250*n_rows),
        grid=false;
        plot_font_size_kwargs...
    )

    plot_idx = 1
    for (network, color) in zip(networks, networks_colors)
        network == "line" && continue
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

        # Calculate mean and std of variance explained for each PC
        max_pcs = maximum(length.(all_variances))
        var_means = zeros(max_pcs)
        var_stds = zeros(max_pcs)
        
        for pc in 1:max_pcs
            pc_vars = [v[pc] for v in all_variances if length(v) >= pc]
            var_means[pc] = mean(pc_vars)
            var_stds[pc] = std(pc_vars)
        end

        # Create subplot
        plot!(
            plt[plot_idx],
            1:max_pcs,
            var_means,
            yerror=var_stds,
            color=color,
            label=nothing,
            title=network,
            xlabel="PC",
            ylabel="Variance explained (%)",
            seriestype=:bar,
            alpha=0.6
        )

        plot_idx += 1
    end

    save_plot(supervisor, plt, "f3_D_intrisic_dim_local_pca")
    display(plt)
end