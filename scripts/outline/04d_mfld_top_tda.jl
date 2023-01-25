using Term.Progress

include("settings.jl")


"""
Load data with single CAN simulations for each network embeded in 10d space, 
run TDA to get topology and local dimensionality estimation.
"""

move_to_datadir(supervisor, "mfld_top")

tda_tresholds = Dict(
    "torus" => 20,
    "sphere" => 18,
    "mobius" => 16,
    "cylinder" => 12,
    "plane" => 25,
)

for network in networks
    network ∉ ("plane",) && continue
    @info "TOPOLOGY analysis: $network"
    filters = Dict{Symbol, Any}(
        :tag => "d10_embeddings",
        :can => network,
    )
    
    do_tda(supervisor, filters, "04_TDA_$(network)"; 
            max_d=2, tresh = tda_tresholds[network],
    )
end

