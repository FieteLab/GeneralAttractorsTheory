using Plots


using GeneralAttractors
using GeneralAttractors.Simulations
using GeneralAttractors.Kernels
import GeneralAttractors.Simulations: plot_trajectory_and_decoded
import GeneralAttractors: animate_simulation_data

include("settings.jl")

"""
Measure error as the distance between the true and decoded trajectory
as given by the CAN's coverspace metric.
"""
function measure_error(can, X, X̄)
    d = metrics[typeof(can.C.M)]

    return [d(X[i, :], X̄[i, :]) for i in 1:size(X, 1)]
end



"""
Make a large number of simulations with random trajectories for 
one CAN and measure error over time averaged across simulations
"""
function run_sims_and_save(network, funky, N_sims, η, duration, still)
    nframes = (Int ∘ round)(duration / dt)
    runs_errors = []
    for i in 1:N_sims
        if i % 10 == 0 || i == 1
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
            σv = 1,
            μv = 0,
            vmax = max_path_int_vel[network],
            still = still,
            x₀ = x₀_traj,
            smoothing_window = 501,
            δ = 10,
            scale=0.5,
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
        error = measure_error(can, trajectory.X, X̄)[still:end]

        # save plots of the trajectory and the decoded trajectory
        if  i % 10 == 0 || i == 1
            trajplot = plot(trajectory; Δ=150)
            display(trajplot)
            save_plot(supervisor, trajplot, "11_path_int_error_$(network)_funky_$(funky)_traj_i_$(i)"; as_svg=false)
        end

        decplt = plot_trajectory_and_decoded(trajectory, X̄) 
        display(decplt)
        (i % 10 == 0 || i == 1) && save_plot(supervisor, decplt, "11_path_int_error_$(network)_funky_$(funky)_traj_i_$(i)"; as_svg=false)

        push!(runs_errors, moving_average(error, 51))
    end




    # ----------------------------------- plot ----------------------------------- #

    E = hcat(runs_errors...)
    store_data(
        supervisor;
        savename = supervisor.datadir/ "path_int_err" / "11_path_int_error_$(network)_funky_$(funky)_η_$(η).npz",
        data = E,
        metadata = Dict(
            "tag" => "path integration error",
            "network" => network,
            "funky" => funky,
            "η" => η,
            "N_sims" => N_sims,
            "duration" => duration,
        )
    )


    # plot error over time
    fig = plot(
        range(0, duration, length=nframes+1),
        mean(E; dims=2),
        ribbon = std(E; dims=2),
        xlabel = "time (ms)",
        ylabel = "on manifold distance (a.u.)",
        title = "Path integration error over time (η = $η - $(network))",
        label = "mean ± std",
        size=(1000, 500);
        plot_font_size_kwargs...
    )

    save_plot(supervisor, fig, "11_path_int_error_$(network)_funky_$(funky)_η_$(η)")
    display(fig)
end


# ---------------------------------------------------------------------------- #
#                                PARAMS & LAUNCH                               #
# ---------------------------------------------------------------------------- #

N_sims = 25

η = 0.0  # noise level
duration = 500

still = 50  # initialization period        
funky = true

# ("plane", "cylinder", "torus")
for network in ("torus",)
    run_sims_and_save(network, funky, N_sims, η, duration, still)
end
