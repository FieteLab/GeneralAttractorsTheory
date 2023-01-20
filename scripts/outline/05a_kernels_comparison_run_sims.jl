using Term.Progress
include("settings.jl")

move_to_datadir(supervisor, "kernels_comparison")

"""
Run a large number of simulations for a single copy CAN with different initial conditions.
Save the data & metadata for further analysis.
"""

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


# # -------------------------------- run & save -------------------------------- #
pbar = ProgressBar()
Progress.with(pbar) do
    job = addjob!(pbar, description = "Simulation", N = N_sims)
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

