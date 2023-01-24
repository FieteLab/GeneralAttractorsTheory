using Term.Progress

include("settings.jl")
import GeneralAttractors.Analysis.ManifoldAnalysis: 
        tda_on_pointcloud, estimate_intrinsic_dimensionality
using Ripserer
import OrderedCollections: OrderedDict
using Term.Tables

"""
Load data with single CAN simulations for each network embeded in 10d space, 
run TDA to get topology and local dimensionality estimation.
"""

move_to_datadir(supervisor, "mfld_top")

tda_params = AnalysisParameters(
    tda_threshold = 30,       # threshold to reduce TDA computation 
    tda_downsample_factor = 10,        # temporal downsampling of data for TDA
    tda_dim_max = 1,        # max feature dimension, starting at 0
)

for network in networks

    # network ∈ ("torus",) && continue

    @info "TOPOLOGY analysis: $network"
    # load data
    filters = Dict{Symbol, Any}(
        :tag => "d10_embeddings",
        :can => network,
        # :dim => 10,
    )
    _, M = ProjectSupervisor.fetch(supervisor; filters...) 
    @assert length(M) == 1 length(M)
    M = M[1]
    @assert size(M, 1) == 10 size(M)

    # run TDA
    tda_model, tda_barcode_plot = tda_on_pointcloud(M, tda_params)
    save_plot(supervisor, tda_barcode_plot, "04_TDA_$(network)");


    # get the number of persistence features based on the largest gap in the persitence
    function get_n_persistence_features(intervals)
        per = persistence.(intervals)
        gaps = per[2:end] .- per[1:end-1]
        largest = argmax(gaps)
        return length(per) - largest
    end
    features = OrderedDict(
        map(d -> d=>get_n_persistence_features(tda_model[d]), 1:tda_params.tda_dim_max+1)...
    )


    Panel(
        Table(Dict(
        "dim" => collect(keys(features)),
        "n_features" => collect(values(features)),
    ));
        title = "$(network)",
        title_style = "blue bold",
        style="bright_yellow",
        fit=false,
        width = 50
    ) |> print

end

