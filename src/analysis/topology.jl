
"""
Function to process (e.g. dimensionality reduction) a point
cloud of data an analyze its topology. 
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

    import GeneralAttractors: load_simulation_history, save_data, load_data, save_model, load_model, savepath
    import GeneralAttractors.Simulations: History
    import GeneralAttractors.Analysis: AnalysisParameters

    export pca_dimensionality_reduction, isomap_dimensionality_reduction, estimate_manifold_topology, estimate_intrinsic_dimensionality

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
            history = load_simulation_history(simulation_name, simulation_name)

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
            save_model(pca_model, simulation_name, "pca_model", :PCA)
            save_data(S_pca_space, simulation_name, "pca_space")
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
        X = load_data(simulation_name, "pca_space")

        # fit
        @info "Performing ISOMAP" size(X) params.n_isomap_dimensions
        iso = ManifoldLearning.fit(
                Isomap,
                X;
                k=params.isomap_k,
                maxoutdim=params.n_isomap_dimensions,
            )
        M = predict(iso, X)
        @info "isomap fitting completed" size(M)

        # make animation of Isomap embedding
        animate_3d_scatter(M, simulation_name*"_iso_embedding"; alpha=.25, title=simulation_name)

        # save
        save_model(iso, simulation_name, "isomap_model", :ISOMAP)
        save_data(M, simulation_name, "isomap_space")
        M
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
    function tda_on_pointcloud(
        simulation_name::String, 
        params::AnalysisParameters,
    )::Vector{PersistenceDiagram}
        # load
        X = load_data(simulation_name, "pca_space")
        tda_on_pointcloud(X, params, simulation_name)
    end

    function tda_on_pointcloud(X::Matrix, params::AnalysisParameters, savename::String)::Vector{PersistenceDiagram}
        # convert M in a vector of tuples for TDA
        n = (Int ∘ round)(size(X, 2)/params.tda_downsample_factor)
        X̄ = [
            Tuple(x) for x in rand(collect(eachcol(X)), n)
        ]

        # fit TDA
        @info "Fitting TDA on X̄" length(X̄) length(X̄[1])
        tda = ripserer(
            X̄; dim_max=params.tda_dim_max, verbose=true, reps=true, threshold=params.tda_threshold
        )

        # plot results
        resplot = plot(
            plot(tda),
            barcode(tda)
        )
        savefig(savepath(savename, "tda_results", "png"))
        save_model(tda, savename, "tda_barcode", :TDA)
        return tda
    end

    function estimate_manifold_topology(simulation_name::String, params::AnalysisParameters=AnalysisParameters())
        # tda = tda_on_pointcloud(simulation_name, params)
        error("Need to analyze the results of TDA on data")
    end

    function estimate_manifold_topology(
                X::Matrix, 
                params::AnalysisParameters=AnalysisParameters(), 
                savename::String="test"
            )
        tda = tda_on_pointcloud(X, params, savename)
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

    Fit PCA to local neighborhoods on the data manifold to estiamte
    local dimensionality. 
    The size of the neighborhood (fraction of datapoints included) and
    the fraction of variance explained requied affectthe results. 
    """
    function estimate_intrinsic_dimensionality(
        simulation_name::String,
        args...;
        kwargs...
    )::Vector{Int}

    # load
    M = load_data(simulation_name, "pca_space")
    estimate_intrinsic_dimensionality(M, args...; kwargs...)
    end
    
    function estimate_intrinsic_dimensionality(
        M::Matrix,
        params::AnalysisParameters=AnalysisParameters(),
        verbose::Bool = true
    )::Vector{Int}

    # build nearest neighbor tree
    nntree = KDTree(M)

    # sample random points on the manifold
    seeds = M[:, rand(1:size(M, 2), params.intrinsic_d_npoints)]

    # get local neighborhoods on the data manifold
    ϵ = params.intrinsic_d_ϵ₀               # search ball radius around each point
    μ, Us, σ = 0.0, nothing, nothing
    while μ < params.intrinsic_d_data_fraction_threshold
        Us = inrange(nntree, seeds, ϵ, true)
        # get coverage
        nᵤ = length.(Us) ./ size(M, 2) * 100
        μ = round(mean(nᵤ); digits=2)
        σ = round(std(nᵤ); digits=2)

        # increment search radius
        ϵ *= 1.5
    end
    verbose && tprintln("$(length(Us)) neighborhoods with $μ{red}±{/red}$σ % of the data each. ($(mean(length.(Us))) data points)")

    # for each neighborhood fit PCA and get the number of PCs
    D = Int[]  # store "dimensionality" of each tangent vector space
    for U in Us
        pca_model = fit(PCA, M[:, U]; pratio=params.intrinsic_d_pratio, );
        push!(D, length(principalvars(pca_model)))
    end
    end

    
end