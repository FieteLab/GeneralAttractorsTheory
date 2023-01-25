
using Term.Progress
include("settings.jl")

move_to_datadir(supervisor, "kernels_comparison")

"""
Similarly to 04a_mfld_top_run_simulations.jl
run a large number of simulations for a single copy 
CAN with different initial conditions and save data and
embeddings to low dimensional spaces.

The only difference is that we only use the torus and 
vary the kernel type.
"""

GENERATE_DATA = false
GENERATE_EMBEDDINGS = true
DO_PLOTS = false

# can params
can_name = "torus"
can_maker = network_makers[can_name]
kernel_params = kernels_parameters_range[can_name]

# number of sims
n_sims_per_kernel = 2500
n_random_kernel_params_draws = 10  # for eacg kernel type
N_sims = n_sims_per_kernel * n_random_kernel_params_draws * length(kernels)
@info "Running $N_sims simulations in total."

# sim params
duration = 70
still = 50
tag = "kernels_comparison"


# -------------------------------- run & save -------------------------------- #
if GENERATE_DATA
    for (name, ktype) in pairs(kernels), rep in 1:n_random_kernel_params_draws
        # make network with random kernel
        params = Dict(k => rand(v) for (k, v) in kernel_params[name])
        kernel = ktype(; params...)
        can = can_maker(:single; k=kernel)

        # run simulations
        for i in 1:n_sims_per_kernel
            generate_or_load(
                supervisor, "simulations";
                fmt = "jld2", 
                name = "history_$(ktype)_rep_$(rep)_$(i)",
                metadata = Dict(
                    :can => can.name,
                    :dt => dt,
                    :duration => duration,
                    :still => still,
                    :τ => τ,
                    :b₀ => b₀,
                    :tag => tag,
                    :n_neurons => collect(can.n),
                    :kernel => string(ktype),
                    :rep => rep,
                )
            ) do 
                h, _ = simulate_constant_traj_random_init(can, duration, dt, still, τ, b₀)
                Dict(h)
            end
            update!(job)
        end
    end
end

# -------------------------------- embeddings -------------------------------- #

if GENERATE_EMBEDDINGS
    for kernel in values(kernels)
        print(hLine(string(kernel); style="red"))
        filters = Dict{Symbol, Any}(
            :tag => "random_initial_conditions",
            :can => can_name,
            :kernel => string(kernel),
        )
        @info "doing kernel $kernel" filters

        metadata = get_entries(supervisor; filters...)
        reps = unique(metadata[:, "rep"]) |> length
        @info "Found $reps reps for kernel `$kernel` simulations."

        for rep in 1:reps
            print(hLine("$rep/$reps"; style="dim blue"))
            # load data for this rep
            filters[:rep] = rep
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
                print(hLine(; style="dim"))
            end
        end
    end
end

# ----------------------------------- plot ----------------------------------- #

if DO_PLOTS
# # M = @time load(; filters...) |> do_pca |> do_isomap

# plots = []
# camera_angles = (45, 50), (20, 20), (70, 70)
# for camera in camera_angles
#     push!(plots, 
#         scatter3d(
#             M[1, :], M[2, :], M[3, :],
#             msa=0, msw=0, alpha=1,
#             markersize = 10.0,
#             legend = false,
#             # title = "Isomap",
#             # xlabel = "M1", ylabel = "M2", zlabel = "M3",
#             color = :black, 
#             camera = camera,
#             size = (800, 800),
#             showaxis = false,
#             axis=nothing,
#         )
#     )
# end
# fig = plot(plots..., layout=(1, 3), size=(2400, 800))
end