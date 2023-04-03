include("settings.jl")
move_to_datadir(supervisor, "mfld_top2")


# done = ("line", "ring", "sphere", "torus", )
todo = ("torus", )

# η = 1.5
# for network in networks

network = "torus"
for η in (0.0, 1.5, 3.0, 5.0)
    (η > 0 && network ∉ ("ring", "torus", "sphere")) && continue
    network ∉ todo && continue

    print(hLine(network * " η:$η"; style="red"))
    
    filters = Dict{Symbol, Any}(
        :tag => "d10_embeddings",
        :can => network,
        :η => η,
        # :funky => false,
    )
    
    max_d = network ∈ ("ring", "line") ? 1 : 2
    do_tda(network, η, supervisor, filters, "f3_C_TDA_$(network)"; 
            max_d=2, tresh = 20,
    )

end

