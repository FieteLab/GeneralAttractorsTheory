
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

    import GeneralAttractors: 
        load_simulation_history,
        save_data,
        load_data,
        save_model,
        load_model,
        savepath,
        checkpath

    import GeneralAttractors.Simulations: History
    import GeneralAttractors.Analysis: AnalysisParameters

    include("analysis_viz.jl")

    export pca_dimensionality_reduction, isomap_dimensionality_reduction, estimate_manifold_topology, estimate_intrinsic_dimensionality

    # ---------------------------------------------------------------------------- #
    #                                      PCA                                     #
    # ---------------------------------------------------------------------------- #
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
        @assert σ[1] > σ[end] string(σ)

        # fit a line through the first and last point
        # consider that σ₁ = f(x₁) - x₁=1
        # we don't start at x=0 so fix the intercept
        m = (σ[end] - σ[1])/(length(σ)-1)
        y₀ = σ[1]
        ŷ(x) = m*x + y₀

        # compute the line between the start and end value
        x = 0:length(σ)-1
        y = ŷ.(x)

        # get the vertical distance between σ and the line
        Δ = σ .- y  

        # get the point with the largest vertical distance
        return argmax(Δ)
    end


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
        )::Nothing
            (checkpath(simulation_name, "pca_space", "npz") && !params.debug) && return

            history = load_simulation_history(simulation_name, simulation_name)

            # flatten out S
            N = size(history.S, 1) * size(history.S, 2)  # tot neurons
            S = reshape(history.S, (N, size(history.S, 3)))[:, 2:end] # N×T | each column an observation

            # get numnber of PCs
            nPC = if !isnothing(params.max_nPC)
                    params.max_nPC
            else
                25
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
            nothing
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
        )::Nothing
            (checkpath(simulation_name, "isomap_space", "npz") && !params.debug) && return

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
            animate_3d_scatter(M, simulation_name, "isomap_projection"; alpha=.25, title=simulation_name)

            # save
            save_model(iso, simulation_name, "isomap_model", :ISOMAP)
            save_data(M, simulation_name, "isomap_space")
            nothing
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

    
    function tda_on_pointcloud(X::Matrix, params::AnalysisParameters, savename::String) #::Tuple{Vector{PersistenceDiagram}, Plot}
        (checkpath(savename, "tda_results", "png") && !params.debug) && return

        # convert M in a vector of tuples for TDA
        n = (Int ∘ round)(size(X, 2)/params.tda_downsample_factor)
        X̄ = [
            Tuple(x) for x in rand(collect(eachcol(X)), n)
        ]

        # fit TDA
        @info "Fitting TDA on X̄" length(X̄) length(X̄[1])
        tda = ripserer(
            X̄; 
            dim_max=params.tda_dim_max, 
            verbose=true, 
            reps=false, 
            threshold=params.tda_threshold
        )

        # plot results
        plt = plot(
            plot(tda),
            barcode(tda),
            size=(1000, 800)
        )
        savefig(savepath(savename, "tda_results", "png"))
        save_model(tda, savename, "tda_barcode", :TDA)
        return tda, plt
    end


    """
    Run TDA and analyze barcode to infer topology.

    TODO: barcode analysis
    """
    function estimate_manifold_topology(
                simulation_name::String,
                params::AnalysisParameters=AnalysisParameters()
            )
        return tda_on_pointcloud(simulation_name, params)
    end

    function estimate_manifold_topology(
                X::Matrix, 
                params::AnalysisParameters=AnalysisParameters(), 
                savename::String="test"
            )
        return tda_on_pointcloud(X, params, savename)
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
        simulation_name::String,
        args...;
        kwargs...
    )

        # load
        M = load_data(simulation_name, "pca_space")
        estimate_intrinsic_dimensionality(M, args...; kwargs...)
    end
    

    function estimate_intrinsic_dimensionality(
            M::Matrix,
            params::AnalysisParameters=AnalysisParameters();
        )
        @info "Estimating intrinsic dimensionality" size(M) params.intrinsic_d_nseeds params.intrinsic_d_neighborhood_size
        
        # build nearest neighbor tree
        nntree = KDTree(M; reorder=false, leafsize=5)

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
            pca_model = fit(PCA, M[:, U]; pratio=0.999999, maxoutdim=size(M, 2))
            d = find_fraction_variance_explained_elbow(principalvars(pca_model))
            push!(D, d)
        end
        return D
    end

    
end