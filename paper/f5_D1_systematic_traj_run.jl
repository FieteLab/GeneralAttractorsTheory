include("settings.jl")
move_to_datadir(supervisor, "TuningCurves")

duration = 12_000 
funky = false
η=0
nframes = (Int ∘ round)(duration / dt)
sampling_n_pts = 31  # adjust this to control grid density
cover_manifold = :default
tag = "PI_systematic_trajectories"


"""
Generate a known fixed trajectory based on the topology
of the variable and neural manifolds defined in the cover space.
"""
function PI_trajectory_maker(can, x₀_traj)
    M, N = typeof(can.C.M), typeof(can.C.N)

    if N == Torus || N == Manifoldℝ² || N == Cylinder
        return make_sp_trajectory(
            can, 
            can.C.M, 
            can.C.N;
            n_pts = sampling_n_pts,  
            n_steps = nframes,
            dt = dt,
            vmax = max_path_int_vel[can.name]
        )
    else
        throw(ArgumentError("make_sp_trajectory is not implemented for $M and $N"))
    end

end

"""
Make a large number of simulations with random trajectories for 
one CAN and save the data
"""
function run_sims_and_save(network, funky, η; cover_manifold=:default)
    savename = "$(network)_systematic_trajectory"
    savename = replace(savename, "." => "_")
    metadata = Dict(
        :can => network,
        :funky => funky,
        :noise => η,
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
        println("running $network")
        # get can
        can, x₀_traj, _ = make_path_int_can(network; funky=funky, random_x0=true, cover_manifold=cover_manifold)

        # get trajectory
        trajectory = PI_trajectory_maker(can, x₀_traj)

        # create trajectory plot
        trajplot = plot(trajectory.X[:, 1], trajectory.X[:, 2])
        display(trajplot)
        
        # get simulation
        simulation = Simulation(can, trajectory; η = η, b₀ = 1.0);
        activate = get_can_initialization_weights(trajectory, can; δ = 0.5)

        # run 
        h, X̄ = @time run_simulation(
            simulation;
            average_over_ms = 0,
            s₀ = 1.0 .* activate,
            discard_first_ms=100,
        )

        # plot decoded
        decplt = plot_trajectory_and_decoded(trajectory, X̄)
        display(decplt)

        # return data to store
        Dict(
            "h" => h, 
            "trajectory" => trajectory, 
            "decoded" => X̄,
            
        )
    end
end
 

for network in networks
    η > 0 && network != "torus" && continue
    funky == true && network ∉ ("torus", "sphere") && continue

    network ∉ ("torus", "cylinder", "plane") && continue
    
    run_sims_and_save(network, funky, η)
end