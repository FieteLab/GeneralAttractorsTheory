include("settings.jl")

"""
Run a large number of simulations for a single copy CAN with different initial conditions.
Save the data & metadata for further analysis.
"""

import GeneralAttractors: r2_to_cy_scaling, cy_to_r2_scaling
move_to_datadir(supervisor, "mfld_top")


GENERATE_DATA = false
GENERATE_EMBEDDINGS = true
GENERATE_DEUBUG_PLOTS = false

# number of sims
n_sims_per_network = 5000
N_sims = n_sims_per_network * length(networks)
GENERATE_DATA && @info "Running $N_sims simulations in total."

# sim params
duration = 25
still = 15
tag = "manifold_topology_data"
kernel_name = :mexican_hat

# -------------------------------- run & save -------------------------------- #

for can_name in networks    
    GENERATE_DATA || break    
    print(hLine(can_name; style="red"))

    # can params
    can_maker = network_makers[can_name]
    kernel_params = Dict(
        k => minimum(v) for (k, v) in kernels_parameters_range[can_name][kernel_name]
    )
    kernel = kernels[kernel_name](; kernel_params...)

    # make network with params
    can = can_maker(:single; k=kernel)

    # run simulations
    for i in 1:n_sims_per_network
        i % 100 == 0 && tprintln("   doing $i/$n_sims_per_network - `$(can_name)`")
        generate_or_load(
            supervisor,
            can.name;
            fmt = "jld2", 
            name = "history_$(can_name)_$(i)",
            metadata = Dict(
                :can => can.name,
                :dt => dt,
                :duration => duration,
                :still => still,
                :τ => τ,
                :b₀ => b₀,
                :tag => tag,
                :n_neurons => collect(can.n),
            ),
            load_existing=false
        ) do 
            h, _ = simulate_constant_traj_random_init(can, duration, dt, still, τ, b₀)
            Dict(h)
        end
    end
end

# ---------------------------- generate embeddings --------------------------- #

GENERATE_EMBEDDINGS && @info "Dimensionality reduction embeddings"
for can_name in networks
    GENERATE_EMBEDDINGS || break
    print(hLine(can_name; style="red"))
    filters = Dict{Symbol, Any}(
        :tag => tag,
        :can => can_name,
    )

    X = load_and_concat_activations(; filters...) 

    # embed in 3 and 10 dimensions
    for (dim, params) in zip(("d3_", "d10_", "d50_"), (dimred_3d_params, dimred_10d_params, dimred_50d_params))
        meta = Dict(
            :can => can_name,
            :dim => parse(Int, dim[2:end-1]),
            :params => Dict(params),
            :tag  => "$(dim)embeddings",
        )

        generate_or_load(
            supervisor, 
            "embeddings";
            name = "$(can_name)_$(dim)",
            fmt = "npz",
            metadata = meta,
            load_existing=false
        ) do
            pca, X̄ = pca_dimensionality_reduction(X, params)
            iso, X̄ = isomap_dimensionality_reduction(X̄, params)
            store_data(supervisor, "embeddings"; name = "$(can_name)_$(dim)embedding", fmt="bson", data=Dict(:iso=>iso, :pca=>pca), metadata = meta)
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

    for can_name in networks
        # plot rand init conditions
        can = network_makers[can_name](:single)
        x₀ = hcat(map(i -> rand(can.C.N), 1:n_sims_per_network)...)
        p1 = scatter(
            eachrow(x₀)...;
            color = :black, ms=25, alpha=.5,
            title = can_name, label=nothing,
        ) 

        # load data and get lattice bump locations
        filters = Dict{Symbol, Any}(
            :tag => tag,
            :can => can_name,
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