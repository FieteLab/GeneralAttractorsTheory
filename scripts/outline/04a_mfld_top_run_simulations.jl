using Term.Progress
include("settings.jl")

"""
Run a large number of simulations for a single copy CAN with different initial conditions.
Save the data & metadata for further analysis.
"""

# can params
can_name = "torus"
can_maker = network_makers[can_name]
kernel_params = kernels_parameters_range[can_name]

# number of sims
n_sims_per_kernel = 500
n_random_kernel_params_draws = 1  # for eacg kernel type
N_sims = n_sims_per_kernel * n_random_kernel_params_draws * length(kernels)
@info "Running $N_sims simulations in total."

# sim params
duration = 125
still = 100
tag = "random_initial_conditions"


# -------------------------------- run & save -------------------------------- #
pbar = ProgressBar()
Progress.with(pbar) do
    job = addjob!(pbar, description = "Simulation", N = N_sims)
    for (name, ktype) in pairs(kernels), rep in 1:n_random_kernel_params_draws

        params = Dict(k => rand(v) for (k, v) in kernel_params[name])
        kernel = ktype(; params...)

        can = can_maker(:single; k=kernel)

        for i in 1:n_sims_per_kernel
            h, _ = simulate_constant_traj_random_init(can, duration, dt, still, τ, b₀)

            store_data(
                supervisor,
                tag,
                can.name;
                data = Dict(h), 
                fmt = "jld2", 
                name = "history_$(ktype)_$(i)",
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
            )
            update!(job)
        end
    end
end

