using Term.Progress
include("settings.jl")

"""
Run a large number of simulations for a single copy CAN with different initial conditions.
Save the data & metadata for further analysis.
"""

move_to_datadir(supervisor, "mfld_top")

GENERATE_DATA = true

# number of sims
n_sims_per_network = 3000
N_sims = n_sims_per_network * length(networks)
GENERATE_DATA && @info "Running $N_sims simulations in total."

# sim params
duration = 120
still = 100
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
            load_existsing=false
        ) do 
            h, _ = simulate_constant_traj_random_init(can, duration, dt, still, τ, b₀)
            Dict(h)
        end
    end
end

# ---------------------------- generate embeddings --------------------------- #

@info "Dimensionality reduction embeddings"
for can_name in networks
    print(hLine(can_name; style="red"))
    filters = Dict{Symbol, Any}(
        :tag => tag,
        :can => can_name,
    )

    X = load_and_concat_activations(; filters...) 

    # embed in 3 and 10 dimensions
    for (dim, params) in zip(("d3_", "d10_"), (dimred_3d_params, dimred_10d_params))

        generate_or_load(
            supervisor, 
            "embeddings";
            name = "$(can_name)_$(dim)",
            fmt = "npz",
            metadata = Dict(
                :can => can_name,
                :dim => parse(Int, dim[2]),
                :params => Dict(params),
                :tag  => "$(dim)embeddings",
            ),
            load_existsing=false
        ) do
            do_isomap(do_pca(X, params), params)
        end
    end
    print(hLine(; style="dim"))
end

