

"""
Similarly to 04_ we run a large number of simulations (on the torus only)
with varying ammount of noise and we repeat the TDA analysis 
to see how the topology changes with noise.
"""


include("settings.jl")
move_to_datadir(supervisor, "mfld_top_noise")


GENERATE_DATA = false
GENERATE_EMBEDDINGS = true


# number of sims
n_sims_per_network = 5000
noise_values = (0.0, 0.5, 1.0)
N_sims = n_sims_per_network * length(noise_values)
GENERATE_DATA && @info "Running $N_sims simulations in total."

# sim params
duration = 25
still = 15
tag = "manifold_topology_vs_noise_data"
kernel_name = :mexican_hat
can_name = "torus"

# -------------------------------- run & save -------------------------------- #
if GENERATE_DATA
    # make network
    can_maker = network_makers[can_name]
    kernel_params = Dict(
        k => mean(v) for (k, v) in kernels_parameters_range[can_name][kernel_name]
    )
    kernel = kernels[kernel_name](; kernel_params...)
    can = can_maker(:single; k=kernel)

    # run simulations
    for η in noise_values
        for i in 1:n_sims_per_network
            i % 100 == 0 && tprintln("   doing $i/$n_sims_per_network - `η = $(η)`")
            generate_or_load(
                supervisor,
                "data";
                fmt = "jld2", 
                name = "history_noie_$(η)_$(i)",
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
end

# ---------------------------- generate embeddings --------------------------- #

GENERATE_EMBEDDINGS && @info "Dimensionality reduction embeddings"
for η in noise_values
    GENERATE_EMBEDDINGS || break
    print(hLine(can_name; style="red"))
    filters = Dict{Symbol, Any}(
        :tag => tag,
        :can => can_name,
        :η => η
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
            load_existing=false
        ) do
            do_isomap(do_pca(X, params), params)
        end
    end
    print(hLine(; style="dim"))
end

# TODO don't embedd all data at once
# TODO do TDA