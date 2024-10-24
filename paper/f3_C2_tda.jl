include("settings.jl")
move_to_datadir(supervisor, "mfld_top4")


network = "klein"
for η in (0.0)  # , 1.5, 3.0, 5.0)
    (η > 0 && network ∉ ("ring", "torus", "sphere")) && continue

    print(hLine(network * " η:$η"; style="red"))
    
    filters = Dict{Symbol, Any}(
        :tag => "d10_embeddings",
        :can => network,
        :η => η,
        :extension => "npz",
        # :funky => false,
    )
    
    max_d = network ∈ ("ring", "line") ? 1 : 2
    do_tda(network, η, supervisor, filters, "f3_C_TDA_$(network)_η_$η"; 
            max_d=2, tresh = 12,  # was 15
    )

end

