using Plots


using GeneralAttractors
using GeneralAttractors.Simulations
using GeneralAttractors.Kernels
import GeneralAttractors.Simulations: plot_trajectory_and_decoded
import GeneralAttractors: animate_simulation_data

include("settings.jl")

"""
Measure error as absolute distance at each time point
"""
function measure_error(X, X̄)
    return sum(abs.(X .- X̄), dims=1) ./ size(X, 1)
end



"""
Make a large number of simulations with random trajectories for 
one CAN and measure error over time averaged across simulations
"""

N_sims = 3

η = 0.0  # noise level
duration = 500
still = 50  # initialization period        
nframes = (Int ∘ round)(duration / dt)

network = "torus"
funky = false

# ------------------------------------ run ----------------------------------- #
runs_errors = []
for i in 1:N_sims
    if i % 10 == 0
        println(
            hLine("run $i / $N_sims - $(network), η=$(η)"; style="bold blue")
        )
    end
    
    # get can
    can, x₀_traj, embedding = make_path_int_can(network; funky=funky, random_x0=true)

    # get trajectory
    trajectory = Trajectory(
        can;
        T = nframes,
        dt = dt,
        σv = 5,
        μv = 0,
        vmax = 0.05,
        still = still,
        x₀ = x₀_traj,
        smoothing_window = 501,
    )

    # get simulation
    simulation = Simulation(can, trajectory; η = 0.0, b₀ = 1.0);
    activate = get_can_initialization_weights(trajectory, can)

    # run 
    h, X̄ = @time run_simulation(
        simulation;
        discard_first_ms = still,
        average_over_ms = 0,
        s₀ = 1.0 .* activate,
    )

    # estimate error
    error = measure_error(trajectory.X, X̄)
    push!(runs_errors, error)
end




# ----------------------------------- plot ----------------------------------- #
# plot error over time
fig = plot(
    range(0, duration, length=nframes),
    mean.(runs_errors),
    ribbon = std.(runs_errors),
    xlabel = "time (ms)",
    ylabel = "abs. error (a.u.)",
    title = "Path integration error over time (η = $η))",
    label = "mean ± std",
    legend = :topleft,
    size = (800, 400),
    dpi = 300,
)

save_plot(supervisor, fig, "11_path_int_error_$(network)_funky_$(funky)_η_$(η)")