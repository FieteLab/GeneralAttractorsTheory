"""
Run a large number of brief simulations with each Single CAN topology
to generate data and embeddings to be used in various figures.

Including with noise. 
"""

include("settings.jl")


move_to_datadir(supervisor, "mfld_top2")




GENERATE_DATA = false
GENERATE_EMBEDDINGS = true
GENERATE_DEUBUG_PLOTS = false

# number of sims
n_sims_per_network = 2500
N_sims = n_sims_per_network * length(networks)
GENERATE_DATA && @info "Running $N_sims simulations in total."

# sim params
duration = 25
still = 15
tag = "manifold_topology_data"

# noise
η = 5.0




# -------------------------------- run & save -------------------------------- #

for network in networks    
    GENERATE_DATA || break    

    (η > 0 && network ∉ ("torus",)) && continue
    print(hLine(network; style="red"))

    # can 
    can = make_single_can(network)
    
    # run simulations
    for i in 1:n_sims_per_network
        i % 100 == 0 && tprintln("   doing $i/$n_sims_per_network - `$(network)`")

        save_name = "$(can.name)_noise_$(η)"
        save_name = replace(save_name, "." => "_")

        name = "history_$(network)_$(i)_noise_$(η)"
        name = replace(name, "." => "_")

        generate_or_load(
            supervisor,
            save_name;
            fmt = "jld2", 
            name = name,
            metadata = Dict(
                :can => can.name,
                :dt => dt,
                :duration => duration,
                :still => still,
                :τ => τ,
                :b₀ => b₀,
                :tag => tag,
                :n_neurons => collect(can.n),
                :η => η,
            ),
            load_existing=false
        ) do 
            h, _ = simulate_constant_traj_random_init(can, duration, dt, still, τ, b₀; η=η)
            Dict(h)
        end
    end
end

# ---------------------------- generate embeddings --------------------------- #

GENERATE_EMBEDDINGS && @info "Dimensionality reduction embeddings"
for network in networks
    (η > 0 && network ∉ ("ring", "torus", "sphere")) && continue
    network != "torus" && continue


    GENERATE_EMBEDDINGS || break
    print(hLine(network * " η:($η)"; style="red"))
    filters = Dict{Symbol, Any}(
        :tag => tag,
        :can => network,
        :η => η,
    )


    X = load_and_concat_activations(; expected_n=n_sims_per_network, filters...) 

    # embed in 3 and 10 dimensions
    for (dim, params) in zip(("d3_", "d10_", "d50_"), (dimred_3d_params, dimred_10d_params, dimred_50d_params))
        meta = Dict(
            :can => network,
            :dim => parse(Int, dim[2:end-1]),
            :params => Dict(params),
            :tag  => "$(dim)embeddings",
            :η => η,
        )

        _η = replace(string(η), "."=>"_")
        name = "$(network)_$(dim)_noise_$(_η)"

        fld = "embeddings_noise_$(_η)"

        generate_or_load(
            supervisor, 
            fld;
            name = name,
            fmt = "npz",
            metadata = meta,
            load_existing=false
        ) do
            pca, X̄ = pca_dimensionality_reduction(X, params)
            iso, X̄ = isomap_dimensionality_reduction(X̄, params)

            # save embeddings in a separate file
            emb_name = "$(network)_$(dim)embedding_noise_$(_η)"
            store_data(supervisor, 
                    "embeddings"; name = emb_name, fmt="bson", data=Dict(:iso=>iso, :pca=>pca), metadata = meta)
            X̄  # save transformed data
        end
    end
    print(hLine(; style="dim"))
end

# ---------------------------- confirmation plots ---------------------------- #

"""
For each network plot the position of initial condition againt
the decoded location of the bump on the manifold
"""

GENERATE_DEUBUG_PLOTS && begin
    @info "Plotting decoded locations"

    for network in networks
        # plot rand init conditions
        can = network_makers[network](:single)
        x₀ = hcat(map(i -> rand(can.C.N), 1:n_sims_per_network)...)
        p1 = scatter(
            eachrow(x₀)...;
            color = :black, ms=25, alpha=.5,
            title = network, label=nothing,
        ) 

        # load data and get lattice bump locations
        filters = Dict{Symbol, Any}(
            :tag => tag,
            :can => network,
        )
        metadata, data = ProjectSupervisor.fetch(supervisor; filters...)
        
        x = hcat(map(
            d -> begin
                    peak = argmax(d["S"][:, 1, end])
                    can.X[:, peak]
                end, data
        )...)

        p2 = scatter(
            eachrow(x)...;
            color = :red, ms=5, alpha=.25, msa=0, msw=0,
            title = "bump locations", label=nothing,
        ) 
        plot(p1, p2) |> display
    end
end