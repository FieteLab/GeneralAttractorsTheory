"""
Trajectories on periodic manifolds plotted in 2
will have "jump" when the trajectory wraps around
a periodic dimension. This removes those jumps
by setting them to NaN so that they don't showup in plots.
"""
function remove_jumps_from_trajectory(X::Matrix)
    Δ = [0, norm.(eachrow(diff(X, dims=1)))...]
    X[Δ .>= 2, :] .*= NaN
    X
end

# ---------------------------------------------------------------------------- #
#                                      viz                                     #
# ---------------------------------------------------------------------------- #


function Plots.plot(traj::Trajectory)
    d = size(traj.X, 2)
    X = remove_jumps_from_trajectory(traj.X)
    if d == 2
        if traj.M isa Mobius
            x, y = X[:, 2], X[:, 1]
        else
            x, y = X[:, 1], X[:, 2]
        end
        plot(
            x,
            y,
            lw = 3,
            color = :black,
            label = nothing,
            title = "trajectory",
            grid = false,
            aspect_ratio = :equal,
        )
    elseif d == 1
        plot(
            X[:, 1],
            lw = 3,
            color = :black,
            label = nothing,
            title = "trajectory",
            grid = false,
        )
    else
        error("Not implemented")
    end
end


function Plots.plot(traj::Trajectory, i::Int; xmin = nothing, xmax = nothing)
    d = size(traj.X, 2)

    if traj.M isa Mobius
        xmin = isnothing(xmin) ? traj.M.xmin : xmin
        xmax = isnothing(xmax) ? traj.M.xmax : xmax
    else    
        xmin = isnothing(xmin) ? minimum(traj.X, dims = 1) : xmin
        xmax = isnothing(xmax) ? maximum(traj.X, dims = 1) : xmax
    end

    t0 = traj.M isa Mobius ? max(1, i-100) : 1
    X = remove_jumps_from_trajectory(traj.X)
    if d == 2
        if traj.M isa Mobius
            k, j = 2, 1
        else
            k, j = 1, 2
        end
        plt = plot(
            X[t0:i, k], X[t0:i, j],
            lw = 4,
            color = :black,
            grid = false,
            aspect_ratio = :equal,
            label = nothing,
            camera = (0.075 * i, 20),
            xlim=[xmin[k], xmax[k]],
            ylim=[xmin[j], xmax[j]],
        )
    elseif d == 1
        plt = plot(
            X[:, 1],
            lw = 4,
            color = :black,
            grid = false,
            label = nothing,
            ylim = [xmin[1], xmax[1]],
        )
    else
        error()
        # plt = plot(
        #     eachcol(X[t0:i, :])...,
        #     lw = 4,
        #     color = :black,
        #     grid = false,
        #     aspect_ratio = :equal,
        #     label = nothing,
        #     _X̄,
        #     ylim = [xmin[2], xmax[2]],
        #     zlim = [xmin[3], xmax[3]],
        #     camera = (0.075 * i, 20),
        # )
    end

    # fain line
    # plot!(plt, eachcol(X)..., lw = 1.5, color = :black, label = nothing, alpha = 0.5)

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

    # th = maximum(s) .* 0.5
    # active = s .>= th
    # inactive = s .< th

    plt = scatter3d(
        eachrow(M)...,
        marker_z = s,
        title = "elapsed: $(round(timems)) ms",
        grid = false,
        msa = 0,
        msw = 0,
        clims = (-maximum(s) * 0.95, maximum(s) * 0.95),
        color=:bwr,
        ms = 6,
        alpha = 0.75,
        label = nothing,
        colorbar = nothing,
    )

    # scatter3d!(
    #     eachrow(M[:, active])...,
    #     marker_z = s[active],
    #     title = "elapsed: $(round(timems)) ms",
    #     grid = false,
    #     msa = 0,
    #     msw = 0,
    #     clims = (0, maximum(s) * 0.95),
    #     ms = 8,
    #     label = nothing,
    #     colorbar = nothing,
    # )

    return plt
end

function custom_sphere_viz(simulation::Simulation, timems, v::Vector, φ; kwargs...)
    can = simulation.can
    s = sum(simulation.S, dims = 2) |> vec

    th = maximum(s) .* 0.5
    active = s .>= th
    inactive = s .< th

    plt = scatter3d(
        eachrow(can.X[:, inactive])...,
        marker_z = s[inactive],
        title = "elapsed: $(round(timems)) ms",
        grid = false,
        msa = 0,
        msw = 0,
        clims = (0, maximum(s) * 0.95),
        ms = 6,
        alpha = 0.1,
        label = nothing,
        colorbar = nothing,
    )
    scatter3d!(
        eachrow(can.X[:, active])...,
        marker_z = s[active],
        title = "elapsed: $(round(timems)) ms",
        grid = false,
        msa = 0,
        msw = 0,
        clims = (0, maximum(s) * 0.95),
        ms = 8,
        label = nothing,
        colorbar = nothing,
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
        end
    else
        error("Simulation plot for d>2 not implemented")
    end


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
        _x = 1:size(simulation.S, 2)
        y = map(i -> maximum(simulation.S[:, i]), _x)
        b = bar(_x, y)

        # visualize distance between decoded and real position
        d = simulation.can.metric(x, X̄[framen, :])
        b2 = bar([0], [d], title = "Dedoded distance", ylim = [0, 3])

        # main figure
        plot(traj, pop_activity, b, b2, size = (1000, 800), layout = (2, 2))
    elseif simulation.can.d == 2
        traj = Plots.plot(tj, framen)

        if simulation.can.name == "mobius"
            i,j = 2,1
        else
            i,j = 1, 2
        end
        
        # plot decoded trajectory
        framen > (100 + 2) && begin
            plot!(
                traj,
                _X̄[100:framen, i],
                _X̄[100:framen, j],
                color = :red,
                label = nothing,
                alpha = 0.6,
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
        plot(traj, pop_activity, size = (1000, 800), layout = (2, 1))
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
function plot_trajectory_and_decoded(trajectory::Trajectory, X̄::Matrix)
    d = size(X̄, 2)
    plt = plot(
        eachcol(trajectory.X)...,
        lw = 5,
        color = :black,
        label = "traj.",
        grid = false,
        aspect_ratio = d > 1 ? :equal : :auto,
        title = "Decoded trajectory",
    )

    scatter!([[x] for x in trajectory.X[1, :]]..., ms = 5, color = :black, label = nothing)

    if d == 2
        plot!(eachcol(X̄)..., lw = 3, color = :red, label = nothing)
    elseif d == 3
        plot3d!(eachcol(X̄)..., lw = 3, color = :red, label = nothing)
    elseif d == 1
        plot!(X̄, lw = 3, color = :red, label = nothing)
    end
    return plt
end





# ---------------------------------------------------------------------------- #
#                                     utils                                    #
# ---------------------------------------------------------------------------- #
reset(simulation::Simulation) = begin
    simulation.Ṡ .*= 0.0
end
reset(Ṡ::Vector) = Ṡ * 0.0
