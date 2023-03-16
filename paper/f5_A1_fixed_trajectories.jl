"""
Run path integration on different "FIXED" trajectories. 
The trajectory differ depending on the topology of the variable manifold.
Different CANs can run the same fixed trajectory. 

Includes running simulations with noise and non-killing vector fields. 
"""

include("settings.jl")
# pyplot()

move_to_datadir(supervisor, "path_int")
tag = "decoding_data"

# TODO add for alternative cover spaces

# ---------------------------------- get CAN --------------------------------- #
duration = 500
still = 150  # initialization period        
funky = false
η = 0.0
nframes = (Int ∘ round)(duration / dt)


# ------------------------ make simulation trajecotry ------------------------ #

function generate_fixed_trajectory(can, network, x₀_traj)
    M, N = typeof(can.C.M), typeof(can.C.N)

    Vs = if M == Manifoldℝ² && N ∈ (Manifoldℝ², Cylinder, Torus)
        v_mag = (cos.(range(0, 2π - .1, length=nframes)) ./ 2 .+ .5)
        vx = sin.(range(0, 1, length=nframes)) .* v_mag ./ 10
        vy = cos.(range(0, 2π - .1, length=nframes)) .* v_mag  ./ 12
        [vx, vy]
        
    elseif M == Cylinder && N == Torus
        x₀_traj = [3.14, -14]
        vx = range(0, 1, length=nframes)  ./ 25
        vy = range(0, .34, length=nframes) .* abs.(1.5 .* sin.(range(0, 6π, length=nframes)) .- .5) ./ 15   
        [vx, vy]

    elseif M == N == Ring
        v_mag = (cos.(range(0, 2π - .1, length=nframes)) ./ 2 .+ .5)
        vx = sin.(range(0, 1, length=nframes)) .* v_mag ./ 10
        [vx]

    elseif M == Line && N == Ring
        v_mag = (cos.(range(0, 5π - .1, length=nframes)) ./ 2 .+ .5)
        vx = sin.(range(-1, 1, length=nframes)) .* v_mag ./ 15
        [vx]

    elseif M == N == Line
        v_mag = (cos.(range(0, 2π - .1, length=nframes)) ./ 2 .+ .5)
        vx = sin.(range(0, .5, length=nframes)) .* v_mag ./ 10
        [vx]

    elseif M ==  N == Mobius
        t  = 1:nframes 
        μ  = abs.(sin.(0.002t) .* .01)

        vx = range(1.5, .5, length = nframes) .* μ
        vy = cos.(0.008t) .* μ ./ 1.25

        σ = 0.1  # scaling
        Δ = 25
        [vx, vy]

    elseif M == Cylinder && N == Mobius
        vx = range(0, 1, length=nframes)  .* .025
        vy = range(0, 1, length=nframes) .* sin.(range(0, 2π, length=nframes)) ./ 10
        [vx, -vy]

    elseif M == N == Sphere
        t  = 1:nframes 
        μ  = abs.(sin.(0.002t) .* .04) .+ .005

        vx = range(.5, 1.0, length = nframes) .* μ
        vy = range(0.5, 1.0, length = nframes) .* μ 
        vz = range(1, 0, length = nframes) .* μ  

        σ = 0.1  # scaling
        Δ = 50
        [σ .* vx, σ .* vy, σ .* vz]


    else
        error("No fixed trajectory for this CAN: $(can)")
    end


    return Trajectory(
        can;
        T = nframes,
        dt = dt,
        vmax = max_path_int_vel[network],
        still = still,
        x₀ = x₀_traj,
        Vs = Vs,
    )
end


function run_network_on_fixed_trajectory(network; savename="", use_x₀=nothing, kwargs...)
    can, x₀_traj, _ = make_path_int_can(network; funky=funky, kwargs...)
    x₀_traj = use_x₀ === nothing ? x₀_traj : use_x₀
    @info "X0" x₀_traj can.C.M can.C.N

    trajectory = generate_fixed_trajectory(can, network, x₀_traj)
    trajplot = plot(trajectory)
    trajplot |> display

    save_plot(supervisor, trajplot, "f5_A1_$(network)_trajectory")
    simulation = Simulation(can, trajectory; η = η, b₀ = b₀);


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
        "η" => η,
    )


    store_data(supervisor, "simulations"; 
            fmt = "jld2", name="$(network)_sim_data_funky_$(funky)_η_$(η)", 
            data = Dict("h" => h, "trajectory"=>trajectory), metadata=meta
    )


    plt = plot_trajectory_and_decoded(trajectory, X̄) 
    plt |> display
    save_plot(supervisor, plt, "f5_A1_fix_traj_PI_$(network)_funky_$(funky)_η_$(η)_decoded_$(savename)")
end

# ---------------------------------------------------------------------------- #
#                                      RUN                                     #
# ---------------------------------------------------------------------------- #

# for network in networks
#     print(hLine(network; style="red"))
#     network != "sphere" && continue

#     η > 0 && network ∉ ("ring", "torus", "sphere") && continue
#     funky && network ∉ ("sphere", "torus") && continue

#     run_network_on_fixed_trajectory(network)
# end


network = "plane"
x0 = [-1.5, 2]
run_network_on_fixed_trajectory(network; 
    # cover_manifold=:line, 
    # savename="alt_cover", 
    # use_x₀=x0 
)
