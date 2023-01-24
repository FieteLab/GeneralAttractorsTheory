using Term.Progress

include("settings.jl")
import GeneralAttractors.Analysis.ManifoldAnalysis: 
        tda_on_pointcloud, estimate_intrinsic_dimensionality


"""
Load data with single CAN simulations for each network embeded in 10d space, 
run TDA to get topology and local dimensionality estimation.
"""

move_to_datadir(supervisor, "mfld_top")

tda_params = AnalysisParameters(
    tda_threshold = 1.25,       # threshold to reduce TDA computation 
    tda_downsample_factor = 7,        # temporal downsampling of data for TDA
    tda_dim_max = 1,        # max feature dimension, starting at 0
    intrinsic_d_nseeds = 200,      # number of seed points for local PCA
    intrinsic_d_neighborhood_size = 200,
)


for network in networks
    # load data
    filters = Dict{Symbol, Any}(
        :tag => "d10_embeddings",
        :can => network,
        :dim => 10,
    )

    _, M = ProjectSupervisor.fetch(supervisor; filters...) 
    @assert length(M) == 1
    M = M[1]
    @assert size(M, 1) == 3

    # run TDA
    tda_model, tda_barcode_plot = tda_on_pointcloud(M, tda_params)
    save_plot(supervisor, tda_barcode_plot, "04_TDA_$(network)");

    # estimate intrinsic dimensionality
    D = estimate_intrinsic_dimensionality(M, tda_params)
    @info "Manifold $(network) has intrinsic dimensionality $(mean(D))±$(std(D))"
end

