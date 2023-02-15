"""
Run a large number of brief simulations with each Single CAN topology
to generate data and embeddings to be used in various figures.

Including with noise. 
"""

include("settings.jl")







GENERATE_DATA = true
GENERATE_EMBEDDINGS = true
GENERATE_DEUBUG_PLOTS = false

# number of sims
n_sims_per_network = 5000
N_sims = n_sims_per_network * length(networks)
GENERATE_DATA && @info "Running $N_sims simulations in total."

# sim params
duration = 50
still = 15
tag = "manifold_topology_data"

η = 0.5

# noise

if η == 0.0
    move_to_datadir(supervisor, "mfld_top")
else
    move_to_datadir(supervisor, "mfld_top_noise")
end

# -------------------------------- run & save -------------------------------- #

for network in networks    
    GENERATE_DATA || break    

    (η > 0 && network ∉ ("ring", "torus", "sphere")) && continue
    print(hLine(network; style="red"))

    # can 
    can = make_single_can(network)
    
    # run simulations
    for i in 1:n_sims_per_network
        i % 100 == 0 && tprintln("   doing $i/$n_sims_per_network - `$(network)`")

        save_name = η == 0.0 ? can.name : "$(can.name)_noise_$(η)"
        save_name = replace(save_name, "." => "_")

        generate_or_load(
            supervisor,
            save_name;
            fmt = "jld2", 
            name = "history_$(network)_$(i)",
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
    GENERATE_EMBEDDINGS || break
    print(hLine(network; style="red"))
    filters = Dict{Symbol, Any}(
        :tag => tag,
        :can => network,
    )

    X = load_and_concat_activations(; filters...) 

    # embed in 3 and 10 dimensions
    for (dim, params) in zip(("d3_", "d10_", "d50_"), (dimred_3d_params, dimred_10d_params, dimred_50d_params))
        meta = Dict(
            :can => network,
            :dim => parse(Int, dim[2:end-1]),
            :params => Dict(params),
            :tag  => "$(dim)embeddings",
        )

        generate_or_load(
            supervisor, 
            "embeddings";
            name = "$(network)_$(dim)",
            fmt = "npz",
            metadata = meta,
            load_existing=false
        ) do
            pca, X̄ = pca_dimensionality_reduction(X, params)
            iso, X̄ = isomap_dimensionality_reduction(X̄, params)
            store_data(supervisor, "embeddings"; name = "$(network)_$(dim)embedding", fmt="bson", data=Dict(:iso=>iso, :pca=>pca), metadata = meta)
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