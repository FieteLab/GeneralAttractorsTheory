using GeneralAttractors
using GeneralAttractors.Simulations
using GeneralAttractors.ProjectSupervisor
import GeneralAttractors.Can: SingleCAN
using Term
using Term.Progress
install_term_stacktrace(hide_frames = true)

using GeneralAttractors.ProjectSupervisor

include("settings.jl")

"""
Run a large number of simulations for a single copy CAN with different initial conditions.
Save the data & metadata for further analysis.
"""


# --------------------------------- functions -------------------------------- #


# ---------------------------------------------------------------------------- #
#                                      RUN                                     #
# ---------------------------------------------------------------------------- #


can = toruscan_single   # type of CAN
N_sims = 5000
dt = 0.5
duration = 125
still = 100
tag = "random_initial_conditions"

@assert can isa SingleCAN "CAN must be a SingleCAN."

pbar = ProgressBar()
Progress.with(pbar) do
    job = addjob!(pbar, description = "Simulation", N = N_sims)
    for i in 1:N_sims
        h, _ = simulate_constant_traj_random_init(can, duration, dt, still, τ, b₀)

        store_data(
            supervisor,
            tag,
            can.name;
            history = (Dict(h), "jld2"),
            metadata = Dict(
                :can => can.name,
                :dt => dt,
                :duration => duration,
                :still => still,
                :τ => τ,
                :b₀ => b₀,
                :x₀ => x₀,
                :tag => tag,
                :n_neurons => collect(can.n),
            )
        )
        update!(job)
    end
end

supervisor