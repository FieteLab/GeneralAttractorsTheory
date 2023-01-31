using Plots


using GeneralAttractors
using GeneralAttractors.Simulations
using GeneralAttractors.Kernels
import GeneralAttractors.Simulations: plot_trajectory_and_decoded
import GeneralAttractors: animate_simulation_data

include("settings.jl")
move_to_datadir(supervisor, "path_int")
tag = "decoding_data"


function do_10()

    # ---------------------------------- get CAN --------------------------------- #

    duration = 500
    still = 50  # initialization period        

    network = "sphere"
    can, x₀_traj, embedding = make_path_int_can(network)


    # ------------------------ make simulation trajecotry ------------------------ #
    nframes = (Int ∘ round)(duration / dt)

    if network == "mobius"
        v_mag = (cos.(range(0, 2π - .1, length=nframes)) ./ 2 .+ .5)
        vx = ones(nframes) / 50 .* v_mag
        vy = cos.(range(0, 15π, length=nframes)) / 150 .* v_mag
        Δ = 25
        vmax, Vs = 0.05, [vx, vy]
    else
        v_mag = (cos.(range(0, 2π - .1, length=nframes)) ./ 2 .+ .5)
        vx = sin.(range(0, 2π, length=nframes)) .* v_mag
        # vy = cos.(range(0, 2π, length=nframes)) .* v_mag
        vy = range(1, -2.5, length=nframes)
        vz = range(-1, 1.0, length=nframes)
        Δ = 50
        vmax, Vs = 0.0075, [vx, vy, vz]
    end


    trajectory = Trajectory(
        can;
        T = nframes,
        dt = dt,
        σv = 5,
        μv = 0,
        vmax = vmax,
        still = still,
        x₀ = x₀_traj,
        smoothing_window = 501,
        Vs = Vs,
    )
    trajplot = plot(trajectory; Δ=Δ)
    save_plot(supervisor, trajplot, "10_path_int_$(network)_traj")
    simulation = Simulation(can, trajectory; η = 0.0, b₀ = 1.0);

    # --------------------------------- simulate --------------------------------- #

    activate = get_can_initialization_weights(trajectory, can; δ=0.25)

    h, X̄ = @time run_simulation(
        simulation;
        discard_first_ms = still,
        average_over_ms = 0,
        s₀ = 1.0 .* activate,
    )


    meta  = Dict(
        "network" => network,
        "dt" => dt,
        "duration" => duration,
        "still" => still,
        "tag" => tag,
    )


    store_data(supervisor, "simulations"; fmt = "jld2", name="$(network)_sim_data", data = Dict("h" => h, "trajectory"=>trajectory), metadata=meta)

    # --------------------------------- visualie --------------------------------- #
    plot_trajectory_and_decoded(trajectory, X̄) |> display
    animate_simulation_data(can, trajectory, h, X̄, embedding, 
            (supervisor.projectdir / "plots" /"path_int_$(network)_sim.gif").path
    )
end

do_10()