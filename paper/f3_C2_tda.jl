include("settings.jl")
move_to_datadir(supervisor, "mfld_top")


# done = ("line", "ring", "sphere", "torus", )
todo = ("torus", "sphere")

η = 0.0
for network in networks
    (η > 0 && network ∉ ("ring", "torus", "sphere")) && continue
    # network ∉ todo && continue

    print(hLine(network; style="red"))
    
    filters = Dict{Symbol, Any}(
        :tag => "d10_embeddings",
        :can => network,
        :η => η
    )
    
    max_d = network ∈ ("ring", "line") ? 1 : 2
    do_tda(network, η, supervisor, filters, "f3_C_TDA_$(network)"; 
            max_d=2, tresh = tda_tresholds[network],
    )

end

