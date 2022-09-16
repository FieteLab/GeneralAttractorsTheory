"""
Manifold analysis using
    - PCA: https://juliastats.org/MultivariateStats.jl/stable/pca/#Principal-Component-Analysis
    - ISOMAP: https://wildart.github.io/ManifoldLearning.jl/stable/isomap/
    - Intrinsic dimensionality estimate based on local PCA approximation of tangent spaces
    - TDA: https://mtsch.github.io/Ripserer.jl/dev/
"""
module Analysis
    using Plots
    import Parameters: @with_kw
    using Term.Progress

    import GeneralAttractors: load_simulation_history, save_data, load_data, save_model, load_model, savepath

    include("analysis_viz.jl")


    export AnalysisParameters
    export pca_dimensionality_reduction, isomap_dimensionality_reduction, estimate_manifold_topology, estimate_intrinsic_dimensionality
    export animate_3d_scatter

    # ---------------------------------------------------------------------------- #
    #                                    PARAMS                                    #
    # ---------------------------------------------------------------------------- #
    @with_kw struct AnalysisParameters
        # --------------------------------- topology --------------------------------- #
        # PCA/isomap
        max_nPC::Union{Nothing, Int}                    = nothing  # max num of PCs
        pca_pratio::Float64                             = .8       # fraction of variance explained
        n_isomap_dimensions::Int                        = 3
        isomap_k::Int                                   = 10
        # TDA
        tda_threshold                                   = 10       # threshold to reduce TDA computation 
        tda_downsample_factor::Int                      = 3        # temporal downsampling of data for TDA
        tda_dim_max::Int                                = 2        # max feature dimension, starting at 0
        # intrinsic dimensionality (local PCA)
        intrinsic_d_npoints::Int                        = 200      # number of seed points for local PCA
        intrinsic_d_bbox_fraction_threshold::Float64    = 1.0  # fraction of data in each neighborhood
        intrinsic_d_ϵ₀::Float64                         = 1.0      # starting ball radius for knn
        intrinsic_d_pratio::Float64                     = .95      # fraction of variance explained for local PCA
        
    end

    include("topology.jl")
    using .ManifoldAnalysis: pca_dimensionality_reduction, isomap_dimensionality_reduction, estimate_manifold_topology, estimate_intrinsic_dimensionality
end