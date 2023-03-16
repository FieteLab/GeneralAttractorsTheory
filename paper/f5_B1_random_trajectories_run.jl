"""
Run N simulations with random trajectories to estimate PI
accuracy.
"""

include("settings.jl")
move_to_datadir(supervisor, "PI")


N_sims = 50
still = 150
duration = 1000 + still
funky = false
η=0.0
nframes = (Int ∘ round)(duration / dt)


cover_manifold = :cylinder

tag = "PI_random_trajectories"


"""
Generate a known fixed trajectory based on the topology
of the variable and neural manifolds defined in the cover space.
"""
function PI_trajectory_maker(can, x₀_traj)
    M, N = typeof(can.C.M), typeof(can.C.N)

    kwargs = Dict(
        :T => nframes,
        :dt => dt,
        :μv => 0,
        :vmax => max_path_int_vel[can.name],
        :still => still,
        :scale=>0.5,
    )

    trajectory = if M == N == Line
        Trajectory(
            can;
            σv = [.5,],  # mobius values
            δ = .25,
            kwargs...
        )
    elseif M == N == Ring
        Trajectory(
            can;
            σv = [.5,],  # mobius values
            δ = 0,
            kwargs...
        )
    elseif M ∈ (Manifoldℝ², Cylinder) && N ∈ (Manifoldℝ², Cylinder, Torus)
        Trajectory(
            can;
            σv = [1.0, 1.0],  # mobius values
            δ = .5,
            x₀ = x₀_traj,
            kwargs...
        )
    elseif M == N == Mobius
        Trajectory(
            can;
            σv = [.75, .2],  # mobius values
            δ = 10,
            kwargs...
        )
    elseif M == N == Sphere
        Trajectory(
            can;
            σv = [.1, .1, .1],  # mobius values
            δ = 0,
            kwargs...
        )
    elseif M == Cylinder && N == Mobius
        Trajectory(
            can;
            σv = [.75, .2],  # mobius values
            δ = 10,
            kwargs...
        )
    else
        error("No trajectory maker for $M and $N")
    end

    return trajectory
end

"""
Make a large number of simulations with random trajectories for 
one CAN and save the data
"""
function run_sims_and_save(network, funky, N_sims, η, still; cover_manifold=:default)
    for i in 1:N_sims
        if i % 10 == 0 || i == 1
            println(
                hLine("run $i / $N_sims - $(network), η=$(η)"; style="bold blue")
            )
        end

        savename = "$(network)_funky_$(funky)_noise_$(η)_$(i)"
        savename = replace(savename, "." => "_")
        metadata = Dict(
            :can => network,
            :funky => funky,
            :noise => η,
            :i => i,
            :tag => tag,
            :cover_manifold => cover_manifold,
        )

        generate_or_load(
            supervisor,
            "PI_$(network)_$(cover_manifold)";
            fmt = "jld2", 
            name = savename,
            metadata = metadata,
            load_existing = false,
        ) do
        
            # get can
            can, x₀_traj, _ = make_path_int_can(network; funky=funky, random_x0=true, cover_manifold=cover_manifold)

            # get trajectory
            trajectory = PI_trajectory_maker(can, x₀_traj)

            # get simulation
            simulation = Simulation(can, trajectory; η = η, b₀ = 1.0);
            activate = get_can_initialization_weights(trajectory, can; δ = 0.5)

            # save plots of the trajectory and the decoded trajectory
            trajplot = plot(trajectory; Δ=75)

            
            # run 
            h, X̄ = @time run_simulation(
                simulation;
                discard_first_ms = still,
                average_over_ms = 0,
                s₀ = 1.0 .* activate,
            )

            # plot decoded
            decplt = plot_trajectory_and_decoded(trajectory, X̄)
            display(decplt)

            # return data to store
            Dict("h" => h, "trajectory"=>trajectory, "decoded" => X̄)
        end
    end
end


for network in networks
    η > 0 && network != "torus" && continue
    funky == true && network ∉ ("torus", "sphere") && continue

    network ∉ ("mobius", ) && continue
    
    run_sims_and_save(network, funky, N_sims, η, still; cover_manifold=cover_manifold)
end