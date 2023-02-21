"""
Run path integration on different "FIXED" trajectories. 
The trajectory differ depending on the topology of the variable manifold.
Different CANs can run the same fixed trajectory. 

Includes running simulations with noise and non-killing vector fields. 
"""

include("settings.jl")

import GeneralAttractors.Simulations: plot_trajectory_and_decoded

move_to_datadir(supervisor, "path_int")
tag = "decoding_data"

# TODO make trajectoriers for the Line and Ring CANs
# TODO add for alternative cover spaces

# ---------------------------------- get CAN --------------------------------- #
duration = 500
still = 150  # initialization period        
funky = false

nframes = (Int ∘ round)(duration / dt)

# ------------------------ make simulation trajecotry ------------------------ #

function generate_fixed_trajectory(can, network, x₀_traj)
    M, N = typeof(can.C.M), typeof(can.C.N)
    trajectory = if M == Manifoldℝ² && N ∈ (Manifoldℝ², Cylinder, Torus)
        v_mag = (cos.(range(0, 2π - .1, length=nframes)) ./ 2 .+ .5)
        vx = sin.(range(0, 1, length=nframes)) .* v_mag ./ 10
        vy = cos.(range(0, 2π - .1, length=nframes)) .* v_mag  ./ 12

        Trajectory(
            can;
            T = nframes,
            dt = dt,
            vmax = max_path_int_vel[network],
            still = still,
            x₀ = x₀_traj,
            Vs = [vx, vy],
        )
    elseif M == N == Ring
        v_mag = (cos.(range(0, 2π - .1, length=nframes)) ./ 2 .+ .5)
        vx = sin.(range(0, 1, length=nframes)) .* v_mag ./ 10

        Trajectory(
            can;
            T = nframes,
            dt = dt,
            vmax = max_path_int_vel[network],
            still = still,
            x₀ = x₀_traj,
            Vs = [vx],
        )

    elseif M == N == Line
        v_mag = (cos.(range(0, 2π - .1, length=nframes)) ./ 2 .+ .5)
        vx = sin.(range(0, .5, length=nframes)) .* v_mag ./ 10

        Trajectory(
            can;
            T = nframes,
            dt = dt,
            vmax = max_path_int_vel[network],
            still = still,
            x₀ = x₀_traj,
            Vs = [vx],
        )

    elseif M ==  N == Mobius
        t  = 1:nframes 
        μ  = abs.(sin.(0.002t) .* .01)

        vx = range(1.5, .5, length = nframes) .* μ
        vy = cos.(0.008t) .* μ ./ 1.25

        σ = 0.1  # scaling
        Δ = 25
        Vs = [vx, vy]

        Trajectory(
            can;
            T = nframes,
            dt = dt,
            vmax = max_path_int_vel[network],
            still = still,
            x₀ = x₀_traj,
            Vs = Vs,
        )

    elseif M == N == Sphere
        t  = 1:nframes 
        μ  = abs.(sin.(0.002t) .* .04) .+ .005

        vx = range(.5, 1.0, length = nframes) .* μ
        vy = range(0.5, 1.0, length = nframes) .* μ 
        vz = range(1, 0, length = nframes) .* μ  

        σ = 0.1  # scaling
        Δ = 50
        Vs = [σ .* vx, σ .* vy, σ .* vz]

        Trajectory(
            can;
            T = nframes,
            dt = dt,
            vmax = max_path_int_vel[network],
            still = still,
            x₀ = x₀_traj,
            Vs = Vs,
        )
    else
        error("No fixed trajectory for this CAN: $(can)")
    end

    return trajectory
end


# TODO see why α not working
function run_network_on_fixed_trajectory(network)
    can, x₀_traj, _ = make_path_int_can(network; funky=funky)

    trajectory = generate_fixed_trajectory(can, network, x₀_traj)
    trajplot = plot(trajectory)

    save_plot(supervisor, trajplot, "f5_fix_traj_PI_$(network)_funky_$(funky)_trajectory")
    simulation = Simulation(can, trajectory; η = 0.0, b₀ = b₀);


    activate = get_can_initialization_weights(trajectory, can)

    h, X̄ = run_simulation(
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
        "M" => string(can.C.M),
        "N" => string(can.C.N),
    )


    store_data(supervisor, "simulations"; 
            fmt = "jld2", name="$(network)_sim_data_funky_$(funky)", 
            data = Dict("h" => h, "trajectory"=>trajectory), metadata=meta
    )


    plot_trajectory_and_decoded(trajectory, X̄) |> display
end


for network in networks
    print(hLine(network; style="red"))
    run_network_on_fixed_trajectory(network)
end