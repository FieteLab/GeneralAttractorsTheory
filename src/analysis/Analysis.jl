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
using Statistics


import ..Simulations: decode_peak_location
include("analysis_viz.jl")
include("utils.jl")

export AnalysisParameters
export pca_dimensionality_reduction,
    isomap_dimensionality_reduction,
    estimate_manifold_topology,
    estimate_intrinsic_dimensionality
export animate_3d_scatter

# ---------------------------------------------------------------------------- #
#                                    PARAMS                                    #
# ---------------------------------------------------------------------------- #
@with_kw struct AnalysisParameters
    debug::Bool = false    # if true all analyses are run, otherwise skip previusly completed steps
    # --------------------------------- topology --------------------------------- #
    # PCA/isomap
    max_nPC::Union{Nothing,Int} = nothing  # max num of PCs
    pca_pratio::Float64 = 0.8       # fraction of variance explained
    n_isomap_dimensions::Int = 3
    isomap_k::Int = 10
    isomap_downsample::Int = 10    # temporal downsample factor for fitting isomap

    # TDA
    tda_threshold = 10       # threshold to reduce TDA computation 
    tda_downsample_factor::Int = 3        # temporal downsampling of data for TDA
    tda_dim_max::Int = 2        # max feature dimension, starting at 0

    # intrinsic dimensionality (local PCA)
    intrinsic_d_nseeds::Int = 200      # number of seed points for local PCA
    intrinsic_d_neighborhood_size::Int = 200      # number of points surrounding each seed to include in PCA          
end

include("topology.jl")
using .ManifoldAnalysis:
    pca_dimensionality_reduction,
    isomap_dimensionality_reduction,
    estimate_manifold_topology,
    estimate_intrinsic_dimensionality
end
