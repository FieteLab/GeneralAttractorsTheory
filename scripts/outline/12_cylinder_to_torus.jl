using Plots


using GeneralAttractors
using GeneralAttractors.Simulations
using GeneralAttractors.Kernels
import GeneralAttractors.Simulations: plot_trajectory_and_decoded
import GeneralAttractors: animate_simulation_data

include("settings.jl")
move_to_datadir(supervisor, "path_int")
tag = "cylinder_to_torus"

function run12()

    # ---------------------------------- get CAN --------------------------------- #
    duration = 2000
    still = 50  # initialization period        

    network = "torus"
    can, x₀_traj, embedding = make_path_int_can(network; funky=false, cover_manifold=:cylinder)

    x₀_traj = [3.14, -14]

    # ------------------------ make simulation trajecotry ------------------------ #
    nframes = (Int ∘ round)(duration / dt)
    v_mag = (cos.(range(0, 2π - .1, length=nframes)) ./ 2 .+ .5)
    # vx = sin.(range(0, 1, length=nframes)) .* v_mag ./ 10
    # vx = collect(1:nframes) ./ nframes .* v_mag
    vx = range(0, 1, length=nframes)  ./ 25
    vy = range(0, .34, length=nframes) .* abs.(1.5 .* sin.(range(0, 6π, length=nframes)) .- .5) ./ 15

    trajectory = Trajectory(
        can;
        T = nframes,
        dt = dt,
        σv = 5,
        μv = 0,
        vmax = max_path_int_vel[network],
        still = still,
        x₀ = x₀_traj,
        smoothing_window = 501,
        Vs = [vx, vy],
    )
    trajplot = plot(trajectory)
    display(trajplot)
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
        "funky" => :false,
    )


    store_data(supervisor, "simulations"; fmt = "jld2", name="cylinder_to_torus", data = Dict("h" => h, "trajectory"=>trajectory), metadata=meta)

    # --------------------------------- visualie --------------------------------- #
    plot_trajectory_and_decoded(trajectory, X̄) |> display
    # animate_simulation_data(can, trajectory, h, X̄, embedding, 
    #         (supervisor.projectdir / "plots" /"$(network)_sim_traj.gif").path
    # )


    # --------------------------------- animation -------------------------------- #
    anim = Animation()

    cpts = by_column(
        cylinder_embedding,
        cylinder_maker(:single;cy_extent=20).X)

    tcan = by_column(
        torus_embedding,
        torus_maker(:single).X)


    for i in 50:20:size(X̄, 1)
        p1 = scatter3d(
            eachrow(cpts)...,
            color = :black, label=nothing,
            alpha=.025, ms=10,
            msa=0, msw=0,
            grid=false,
            showaxis = false,
            axis=nothing,
            xlim=[-.8, .8], ylim=[-.8, .8], zlim=[-20, 20]
        )

        p2 = scatter3d(
            eachrow(tcan)...,
            color = :black, label=nothing,
            alpha=.025, ms=10,
            msa=0, msw=0,
            grid=false,
            showaxis = false,
            axis=nothing,
            xlim=[-1.25, 1.25], ylim=[-1.25, 1.25], zlim=[-.8, .8]
        )

        cx = cylinder_embedding(trajectory.X[i, :])
        scatter3d!( p1,
            [cx[1]], [cx[2]], [cx[3]],
            color = :red, ms=10, label=nothing,
            msa=1, msw=1, msc=:white,
        )

        tx = torus_embedding(X̄[i, :])
        scatter3d!( p2,
            [tx[1]], [tx[2]], [tx[3]],
            color = :red, ms=10, label=nothing,
            msa=1, msw=1, msc=:white,
        )

        plot!(p1, p2, size=(800, 800), legend=:none)
        frame(anim)
    end
    gif(anim, "test.gif", fps = 30)

end

run12()