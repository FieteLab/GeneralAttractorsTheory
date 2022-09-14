"""
Manifold analysis using
    - PCA: https://juliastats.org/MultivariateStats.jl/stable/pca/#Principal-Component-Analysis
    - ISOMAP: https://wildart.github.io/ManifoldLearning.jl/stable/isomap/
    - Intrinsic dimensionality estimate based on local PCA approximation of tangent spaces
    - TDA
"""
module Analysis
    using ManifoldLearning
    using NearestNeighbors
    import NearestNeighbors: NNTree
    using MultivariateStats
    import Parameters: @with_kw
    using Plots
    using Term.Progress

    include("analysis_viz.jl")

    import GeneralAttractors: load_simulation_history, save_data, load_data, save_model, load_model, savepath
    import ..Simulations: History

    export pca_dimensionality_reduction, isomap_dimensionality_reduction, AnalysisParameters
    export animate_3d_scatter

    # ---------------------------------------------------------------------------- #
    #                                    PARAMS                                    #
    # ---------------------------------------------------------------------------- #
    @with_kw struct AnalysisParameters
        max_nPC::Union{Nothing, Int}    = nothing  # max num of PCs
        pca_pratio::Float64             = .8       # fraction of variance explained
        n_isomap_dimensions::Int        = 3
        isomap_k::Int                   = 10
        nntype                          = KDTree  # used for isomap and local PCA | <:NNTree
    end


    # ---------------------------------------------------------------------------- #
    #                                      PCA                                     #
    # ---------------------------------------------------------------------------- #
    fraction_variance_explained(M::PCA) = principalvars(M) ./ tvar(M) * 100


    """
        pca_dimensionality_reduction(
            simulation_name::String, 
            params::AnalysisParameters=AnalysisParameters(),
        )

        Reduce dimensionality of simulation data through PCA.
        It assumes a file "./data/simulation_name.bson" exists
        and loads its content to perform PCA.
        Then saves "./data/simulation_name_pca.bson" with the
        PCA projected data.
    """
    function pca_dimensionality_reduction(
            simulation_name::String, 
            params::AnalysisParameters=AnalysisParameters(),
        )   
            history = load_simulation_history(simulation_name)

            # flatten out S
            N = size(history.S, 1) * size(history.S, 2)  # tot neurons
            S = reshape(history.S, (N, size(history.S, 3)))[:, 2:end] # N×T | each column an observation

            # get numnber of PCs
            nPC = if !isnothing(params.max_nPC)
                    params.max_nPC
            else
                (Int ∘ round)(max(N/100, 50))
            end

            # reduce dimensionality with PCA
            @info "PCA dimensionality reduction" size(S) params.max_nPC params.pca_pratio nPC
            pca_model = fit(PCA, S; pratio=params.pca_pratio, maxoutdim=nPC);
            S_pca_space = predict(pca_model, S)
            @info "pca fitting completed $(length(principalvars(pca_model))) PCs"  size(S_pca_space)

            # plot fraction of variance explained
            plot(
                cumsum(fraction_variance_explained(pca_model)),
                legend=nothing, ylabel="cum frac var", xlabel="PC",
                lw=4, color=:black, title="Fraction of variance explained"
            ) |> display

            # save results to file
            save_model(pca_model, simulation_name*"_pca_model", :PCA)
            save_data(S_pca_space, simulation_name*"_pca")
    end


    # ---------------------------------------------------------------------------- #
    #                                    ISOMAP                                    #
    # ---------------------------------------------------------------------------- #
    """
    isomap_dimensionality_reduction(
            simulation_name::String, 
            params::AnalysisParameters=AnalysisParameters(),
        )  

        Performs ISOMAP dimensionality reduction on data that previously
        underwent PCA dimensionality reduction
    """
    function isomap_dimensionality_reduction(
        simulation_name::String, 
        params::AnalysisParameters=AnalysisParameters(),
    )   
        # load
        X = load_data(simulation_name*"_pca")

        # fit
        @info "Performing ISOMAP" size(X) params.n_isomap_dimensions
        iso = ManifoldLearning.fit(
                Isomap,
                X;
                k=params.isomap_k,
                maxoutdim=params.n_isomap_dimensions,
                # nntype=params.nntype
            )
        M = predict(iso, X)
        @info "isomap fitting completed" size(M)

        # make animation of Isomap embedding
        animate_3d_scatter(M, simulation_name*"_iso_embedding"; alpha=.25, title=simulation_name)

        # save
        save_model(iso, simulation_name*"_iso_model", :ISOMAP)
        save_data(M, simulation_name*"_iso")
        M
    end
end