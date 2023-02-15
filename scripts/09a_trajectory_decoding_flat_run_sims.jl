using Plots


using GeneralAttractors
using GeneralAttractors.Simulations
using GeneralAttractors.Kernels
import GeneralAttractors.Simulations: plot_trajectory_and_decoded
import GeneralAttractors: animate_simulation_data

include("settings.jl")
move_to_datadir(supervisor, "path_int")
tag = "decoding_data"

function run09()

    # ---------------------------------- get CAN --------------------------------- #
    duration = 500
    still = 50  # initialization period        

    network = "torus"
    funky = true
    can, x₀_traj, embedding = make_path_int_can(network; funky=funky)

    # ------------------------ make simulation trajecotry ------------------------ #
    nframes = (Int ∘ round)(duration / dt)
    v_mag = (cos.(range(0, 2π - .1, length=nframes)) ./ 2 .+ .5)
    vx = sin.(range(0, 1, length=nframes)) .* v_mag ./ 10
    # vx = collect(1:nframes) ./ nframes .* v_mag
    vy = cos.(range(0, 2π - .1, length=nframes)) .* v_mag  ./ 12

    trajectory = Trajectory(
        can;
        T = nframes,
        dt = dt,
        σv = 5,
        μv = 0,
        vmax = max_path_int_vel[network],
        still = still,
        x₀ = x₀_traj,
        Vs = [vx, vy],
    )
    trajplot = plot(trajectory)
    save_plot(supervisor, trajplot, "09_path_int_traj")
    simulation = Simulation(can, trajectory; η = 0.0, b₀ = 1.0);


    # --------------------------------- simulate --------------------------------- #

    activate = get_can_initialization_weights(trajectory, can)

    h, X̄ = @time run_simulation(
        simulation;
        discard_first_ms = still,
        average_over_ms = 0,
        s₀ = 1.0 .* activate,
    )


    meta = Dict(
        "network" => network,
        "dt" => dt,
        "duration" => duration,
        "still" => still,
        "tag" => tag,
        "funky" => funky,
    )


    store_data(supervisor, "simulations"; fmt = "jld2", name="$(network)_sim_data_funky_$(funky)", data = Dict("h" => h, "trajectory"=>trajectory), metadata=meta)

    # --------------------------------- visualie --------------------------------- #
    plot_trajectory_and_decoded(trajectory, X̄) |> display
    # animate_simulation_data(can, trajectory, h, X̄, embedding, 
    #         (supervisor.projectdir / "plots" /"$(network)_sim_traj.gif").path
    # )
end

run09()