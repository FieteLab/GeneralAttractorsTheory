include("settings.jl")



for network in networks
    @info "TOPOLOGY analysis: $network"
    filters = Dict{Symbol, Any}(
        :tag => "d10_embeddings",
        :can => network,
    )
    
    do_tda(supervisor, filters, "f3_C_TDA_$(network)"; 
            max_d=2, tresh = tda_tresholds[network],
    )
end

