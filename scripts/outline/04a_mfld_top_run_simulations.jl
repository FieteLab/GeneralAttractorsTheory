using Term.Progress
include("settings.jl")

"""
Run a large number of simulations for a single copy CAN with different initial conditions.
Save the data & metadata for further analysis.
"""

move_to_datadir(supervisor, "mfld_top")


# number of sims
n_sims_per_network = 2500
N_sims = n_sims_per_network * length(networks)
@info "Running $N_sims simulations in total."

# sim params
duration = 120
still = 100
tag = "manifold_topology_data"
kernel_name = :local_global

# -------------------------------- run & save -------------------------------- #
pbar = ProgressBar()
Progress.with(pbar) do
    job = addjob!(pbar, description = "Simulation", N = N_sims)
    for can_name in networks        
        # can params
        can_maker = network_makers[can_name]
        kernel_params = Dict(
            k => mean(v) for (k, v) in kernels_parameters_range[can_name][kernel_name]
        )
        kernel = LocalGlobalKernel(; kernel_params...)

        # make network with params
        can = can_maker(:single; k=kernel)

        # run simulations
        for i in 1:n_sims_per_network
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
                )
            ) do 
                h, _ = simulate_constant_traj_random_init(can, duration, dt, still, τ, b₀)
                Dict(h)
            end
            update!(job)
        end
    end
end

# ---------------------------- generate embeddings --------------------------- #

@info "Dimensionality reduction embeddings"
for can_name in networks
    print(hLine(can_name; style="red"))
    filters = Dict{Symbol, Any}(
        :tag => "manifold_topology_data",
        :can => can_name,
    )

    X = load(; filters...) 

    # embed in 3 and 10 dimensions
    for (dim, params) in zip(("d3_", "d10_"), (params_3d, params_10d))

        generate_or_load(
            supervisor, 
            "random_initial_conditions",
            "embeddings";
            name = "$(kernel)_$(dim)_rep_$(rep)",
            fmt = "npz",
            metadata = Dict(
                :can => can_name,
                :kernel => string(kernel),
                :rep => rep,
                :params => Dict(params),
                :tag  => "ric_embeddings",
            )
        ) do
            M = do_isomap(do_pca(X, params), params)
        end
    end
    print(hLine(; style="dim"))
end

