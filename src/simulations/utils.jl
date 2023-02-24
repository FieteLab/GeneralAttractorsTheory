"""
Trajectories on periodic manifolds plotted in 2
will have "jump" when the trajectory wraps around
a periodic dimension. This removes those jumps
by setting them to NaN so that they don't showup in plots.
"""
function remove_jumps_from_trajectory(X::Matrix)
    Δ = [0, norm.(eachrow(diff(X, dims = 1)))...]
    X[Δ.>=2, :] .*= NaN
    X
end

# ---------------------------------------------------------------------------- #
#                                      viz                                     #
# ---------------------------------------------------------------------------- #


function Plots.plot(traj::Trajectory; Δ=50)
    X = remove_jumps_from_trajectory(traj.X)
    X̄ = remove_jumps_from_trajectory(traj.X̄)

    if size(X, 2) == 1
        plt = plot(X[traj.still:end], lw=5, color=:black, label="input trajectory")
        plot!(X̄, lw=5, color=:red, label="output trajectory")
        return plt
    end


    zlim = size(traj.X, 2) > 2 ? [traj.M.xmin[3], traj.M.xmax[3]] : [-1, 1]
 
    p1 = plot(
        eachcol(X)...,
        lw = 3,
        color = :black,
        label = nothing,
        title = "trajectory on: $(traj.M.name)",
        grid = false,
        aspect_ratio = :equal,
        xlim = [traj.M.xmin[1], traj.M.xmax[1]],
        ylim = [traj.M.xmin[2], traj.M.xmax[2]],
        # zlim=zlim,
        right_margin = 12Plots.mm,
        left_margin = 12Plots.mm,
        top_margin = 12Plots.mm,
        bottom_margin = 12Plots.mm,
        dpi=400,
    )

    scatter!(
        eachcol(X[1:Δ:end, :])...,
        marker_z = norm.(eachrow(traj.V))[1:Δ:end],
        # marker_z = 1:size(X[1:Δ:end, :], 1),
        label = nothing,
        msw = .5, msa = .5,
        colorbar_title = "speed",
    )

    p2 = plot(
        eachcol(X̄)...,
        lw = 3,
        color = :red,
        label = nothing,
        title = "trjectory on: $(traj.N.name)",
        grid = false,
        aspect_ratio = :equal,
        xlim = [traj.N.xmin[1], traj.N.xmax[1]],
        ylim = [traj.N.xmin[2], traj.N.xmax[2]],
    )

    plot(p1, p2; size=(1000, 1000))
end


function Plots.plot(traj::Trajectory, i::Int; kwargs...)
    d = size(traj.X, 2)

    lims = d == 2 ? 
        Dict(
        :xlim => [traj.M.xmin[1], traj.M.xmax[1]],
        :ylim => [traj.M.xmin[2], traj.M.xmax[2]],
        ) : Dict(
            :xlim => (-1.1, 1.1),
            :ylim => (-1.1, 1.1),
            :zlim => (-1.1, 1.1),
        )

         

    X = remove_jumps_from_trajectory(traj.X)
    plt = plot(
        eachcol(X[1:i, :])...,
        lw = 6,
        color = :black,
        grid = false,
        aspect_ratio = :equal,
        label = "input trajectory",
        camera = (0.01 * i, 20);
        lims...
    )

    plt
end

# -------------------------------- simulation -------------------------------- #

function simulation_frame_1dcan(simulation::Simulation, timems, v::Vector; kwargs...)
    plt = plot(;
        title = "elapsed: $(round(timems)) ms",
        ylim = [min(-0.1, minimum(simulation.S)), max(0.4, maximum(simulation.S))],
        grid = false,
    )

    plot!(simulation.can.X[1, :], simulation.S[:, 1], lw = 3, label = nothing)
    plot!(simulation.can.X[1, :], simulation.S[:, 2], lw = 3, label = nothing)
    plt
end


function simulation_frame_2dcan(
    simulation::Simulation,
    timems,
    v::Vector,
    ::Nothing;
    kwargs...,
)
    can = simulation.can
    # plt = plot()
    # for i in 1:size(simulation.S, 2)
    #     s = sum(reshape(simulation.S[:, i], can.n), dims=2)
    #     plot!(s, lw=4, label="Pop: $i", alpha=.5)
    # end
    s̄ = sum(simulation.S, dims = 2) |> vec

    plt = plot(;
        title = "elapsed: $(round(timems)) ms",
        clims = (min(0, minimum(s̄)), max(maximum(s̄) / 2, 0.1)),
        aspect_ratio = :equal,
        grid = false,
        size = simulation.can.n .* 10,
    )

    h = maximum(can.X) / 2.5
    v̂ = v ./ norm(v) .* h
    x̄ = range(minimum(can.X[1, :]), maximum(can.X[1, :]), length = can.n[1])
    ȳ = range(minimum(can.X[2, :]), maximum(can.X[2, :]), length = can.n[2])

    contourf!(x̄, ȳ, reshape(s̄, can.n)', levels = 3)
    plt
end


function simulation_frame_2dcan(
    simulation::Simulation,
    timems,
    v::Vector,
    φ::Function;
    kwargs...,
)
    can = simulation.can
    s = sum(simulation.S, dims = 2) |> vec
    M = by_column(φ, can.X)

    plt = scatter3d(
        eachrow(M)...,
        marker_z = s,
        title = "elapsed: $(round(timems)) ms",
        grid = false,
        msa = 0,
        msw = 0,
        clims = (-maximum(s) * 0.95, maximum(s) * 0.95),
        color = :bwr,
        ms = 6,
        alpha = 0.8,
        label = nothing,
        colorbar = nothing,
        camera = (90, 0),
    )

    return plt
end

function custom_sphere_viz(simulation::Simulation, timems, v::Vector, φ; kwargs...)
    can = simulation.can
    s = sum(simulation.S, dims = 2) |> vec


    plt = scatter3d(
        eachrow(can.X)...,
        marker_z = s,
        title = "elapsed: $(round(timems)) ms",
        grid = false,
        msa = 0.4,
        msw = 0.5,
        clims = (-maximum(s) * 0.95, maximum(s) * 0.95),
        color = :bwr,
        ms = 3,
        alpha = 0.9,
        label = nothing,
        colorbar = nothing,
        xlim = [-1.1, 1.1],
        ylim = [-1.1, 1.1],
        zlim = [-1.1, 1.1],
    )

    return plt
end



function Plots.plot(
    simulation::Simulation,
    timems,
    framen,
    x::Vector,   # trajectory's position
    v::Vector,
    X̄,
    φ::Union{Nothing,Function};           # decoded position
    kwargs...,
)
    # plot populatiuon actvitiy
    if simulation.can.d == 1
        pop_activity = simulation_frame_1dcan(simulation, timems, v; kwargs...)
    elseif simulation.can.d == 2
        if simulation.can.name == "sphere"
            pop_activity = custom_sphere_viz(simulation, timems, v, φ; kwargs...)
        else
            pop_activity = simulation_frame_2dcan(simulation, timems, v, φ; kwargs...)

            isa(simulation.can.C.M, Mobius) || begin
                traj_on_mf = remove_jumps_from_trajectory(simulation.trajectory.X̄)
                plot!(
                    pop_activity,
                    traj_on_mf[:, 1], traj_on_mf[:, 2],
                    lw = 3,
                    color = :red,
                    label = nothing,
                )
            end
        end
    else
        error("Simulation plot for d>2 not implemented")
    end

    # plot vectori field over population activity
    plot_can_vector_fields!(pop_activity, simulation.can, v, x, X̄[framen, :])


    # plot trajectory
    tj = simulation.trajectory
    # ϕ(x) = moving_average(x, 31)
    _X̄ = remove_jumps_from_trajectory(X̄)
    if simulation.can.name == "sphere"
        traj = Plots.plot(tj, framen;)
        scatter3d!(traj, [x[1]], [x[2]], [x[3]], ms = 5, color = :black, label = "actual")

        # plot decoded trajectory
        plot3d!(
            traj,
            _X̄[1:framen, 1],
            _X̄[1:framen, 2],
            _X̄[1:framen, 3],
            color = :red,
            label = nothing,
            alpha = 1,
            lw = 2,
        )
        scatter3d!(
            traj,
            [_X̄[framen, 1]],
            [_X̄[framen, 2]],
            [_X̄[framen, 3]],
            ms = 7,
            color = :red,
            label = nothing,
        )

        # visualize activation of each individual population
        # _x = 1:size(simulation.S, 2)
        # y = map(i -> maximum(simulation.S[:, i]), _x)
        # b = bar(_x, y)

        # # visualize distance between decoded and real position
        # d = simulation.can.metric(x, X̄[framen, :])
        # b2 = bar([0], [d], title = "Dedoded distance", ylim = [0, 3])

        # # main figure
        # plot(traj, pop_activity, b, b2, size = (1000, 800), layout = (2, 2))
        plot(traj, pop_activity)
    elseif simulation.can.d == 2
        traj = Plots.plot(tj, framen)

        # i, j = simulation.can.name == "mobius" ? (2, 1) : (1, 2)
        i, j = 1, 2

        # plot decoded trajectory
        framen > (50 + 2) && begin
            plot!(
                traj,
                _X̄[50:framen, i],
                _X̄[50:framen, j],
                color = :red,
                label = nothing,
                alpha = 0.6,
                lw = 4,
            )
        end
        scatter!(
            traj,
            [_X̄[framen, i]],
            [_X̄[framen, j]],
            ms = 7,
            color = :red,
            label = nothing,
        )

        # main figure
        plot(traj, pop_activity, size = (1000, 800), layout = (1, 2))
    elseif simulation.can.d == 1
        traj = Plots.plot(tj, framen)

        # plot decoded trajectory
        framen > (100 + 2) &&
            plot!(traj, _X̄[1:framen, 1], color = :red, label = nothing, alpha = 0.6)
        scatter!(traj, [framen], [_X̄[framen, 1]], ms = 7, color = :red, label = nothing)

        # main figure
        plot(traj, pop_activity, size = (1000, 800), layout = (2, 1), ylim = [0, 2π])
    else
        traj = Plots.plot(tj, framen)

        # plot decoded trajectory
        framen > (100 + 2) && plot!(
            traj,
            _X̄[100:framen, 1],
            _X̄[100:framen, 2],
            color = :red,
            label = nothing,
            alpha = 0.6,
        )
        scatter!(
            traj,
            [_X̄[framen, 1]],
            [_X̄[framen, 2]],
            ms = 7,
            color = :red,
            label = nothing,
        )

        # main figure
        plot(traj, pop_activity, size = (1000, 800), layout = (1, 2))
    end
end


# ---------------------------------------------------------------------------- #
#                               END OF SIMULATION                              #
# ---------------------------------------------------------------------------- #
function plot_trajectory_and_decoded(trajectory::Trajectory, X̄::Matrix; kwargs...)
    X = remove_jumps_from_trajectory(trajectory.X)
    X̄ = remove_jumps_from_trajectory(X̄)

    if size(X, 2) == 1
        plt = plot(X, color=:black, label="traj.")
        plot!(X̄, color=:red, label="decoded")
        return plt
    end

    d = size(X̄, 2)
    plt = plot(
        eachcol(X)...,
        lw = 5,
        color = :black,
        label = "traj.",
        grid = false,
        aspect_ratio = d > 1 ? :equal : :auto,
        title = "Decoded trajectory",
        size=(1000, 1000);
        kwargs...
    )

    # scatter!([[x] for x in X[1, :]]..., ms = 5, alpha=.2, color = :black, label = nothing)
    scatter!([[x] for x in X[end, :]]..., ms = 5, alpha=.8, color = :white, msc=:black, msw=2, label = nothing)

    if d == 2
        plot!(eachcol(X̄)..., lw = 3, color = :red, label = nothing)
        # scatter!([[x] for x in X̄[1, :]]..., ms = 5, alpha=.8, color = :red, label = nothing)
        scatter!([[x] for x in X̄[end, :]]..., ms = 5, alpha=.8, color = :white, msc=:red, msw=2, label = nothing)
    elseif d == 3
        plot3d!(
            eachcol(X̄)...,
            lw = 3,
            color = :red,
            label = nothing,
            xlim = [-1.1,1.1],
            ylim = [-1.1, 1.1],
            zlim = [-1.1, 1.1],
            grid = false,
            camera = (190, 30)
        )
    elseif d == 1
        plot!(X̄, lw = 3, color = :red, label = nothing)
    end
    return plt
end


"""
Plot a trajectory's on-manifold trace and the neural data's one.
"""
function plot_on_mfld_trajectory_and_history(can, trajectory::Trajectory, h::History)
    i, j = trajectory.M isa Mobius ? (2, 1) : (1, 2)
    X̄ = remove_jumps_from_trajectory(trajectory.X̄)

    p = plot(
            X̄[:, i], X̄[:, j],
            lw = 3,
            color = :black,
            label = nothing,
            grid = false,
            aspect_ratio = :equal,
            title = "on manifold trajectory"
        )

    S = h.S
    n, _, m = size(S)
    S = real.(reshape(sum(S, dims = 2), (n, m)))
    peak_location =
        Matrix(hcat(map(st -> decode_peak_location(st, can), eachcol(S))...)') |>
        remove_jumps_from_trajectory

    plot!(
        p,
        peak_location[:, i],
        peak_location[:, j],
        lw = 3,
        color = :red,
        label = nothing,
    )
    p
end





# ---------------------------------------------------------------------------- #
#                                     utils                                    #
# ---------------------------------------------------------------------------- #
reset(simulation::Simulation) = begin
    simulation.Ṡ .*= 0.0
end
reset(Ṡ::Vector) = Ṡ * 0.0
