using GeneralAttractors
using GeneralAttractors.Simulations
using GeneralAttractors.ProjectSupervisor
import GeneralAttractors.Can: SingleCAN
using Term
using Term.Progress
install_term_stacktrace(hide_frames = true)



"""
Run a large number of simulations for a single copy CAN with different initial conditions.
Save the data & metadata for further analysis.
"""

supervisor = Supervisor("GeneralAttractorsTheory")
set_datadir(supervisor, "/Users/federicoclaudi/Dropbox (Personal)/Postdoc/GeneralAttractors/data")

# --------------------------------- functions -------------------------------- #
function random_init(can)
    x₀ = can.C.N |> rand
    d = map(i -> can.metric(x₀, can.X[:, i]), 1:size(can.X, 2))
    activate = zeros(length(d))
    activate[d.<0.5] .= 1
    return x₀, activate
end

function setup(can, duration, dt, still, τ, b₀, x₀)
    # initialize trajectory and simulation
    nframes = (Int ∘ round)(duration / dt)
    trajectory = ConstantTrajectory(
        can;
        T = nframes,
        still = still,
    )

    return Simulation(can, trajectory; η = 0.0, b₀ = b₀, τ = τ)
end


# ---------------------------------------------------------------------------- #
#                                      RUN                                     #
# ---------------------------------------------------------------------------- #


can = toruscan_single   # type of CAN
N_sims = 5000
dt = 0.5
duration = 125
still = 100
τ = 5.0
b₀ = 0.5
tag = "random_initial_conditions"

@assert can isa SingleCAN "CAN must be a SingleCAN."

pbar = ProgressBar()
Progress.with(pbar) do
    job = addjob!(pbar, description = "Simulation", N = N_sims)
    for i in 1:N_sims
        x₀, activate = random_init(can)
        sim = setup(can, duration, dt, still, τ, b₀, x₀)
        h, X = run_simulation(    
            sim;
            discard_first_ms = still,
            average_over_ms = 1,
            s₀ = 1.0 .* activate,
        );

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