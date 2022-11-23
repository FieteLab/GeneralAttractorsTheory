# ---------------------------------------------------------------------------- #
#                                      viz                                     #
# ---------------------------------------------------------------------------- #


function Plots.plot(traj::Trajectory)
    d = size(traj.X, 2)
    d != 2 && error("not implemented")

    plot(
        traj.X[:, 1],
        traj.X[:, 2],
        lw = 3,
        color = :black,
        label = nothing,
        title = "trajectory",
        grid = false,
        aspect_ratio = :equal,
    )
end


function Plots.plot(traj::Trajectory, i::Int; xmin = nothing, xmax = nothing)
    d = size(traj.X, 2)

    xmin = isnothing(xmin) ? minimum(traj.X, dims = 1) : xmin
    xmax = isnothing(xmax) ? maximum(traj.X, dims = 1) : xmax

    if d == 2
        plt = plot(
            eachcol(traj.X[1:i, :])...,
            lw = 4,
            color = :black,
            grid = false,
            aspect_ratio = :equal,
            label = "trajectory",
            camera = (0.075 * i, 20),
        )
    else
        plt = plot(
            eachcol(traj.X[1:i, :])...,
            lw = 4,
            color = :black,
            grid = false,
            aspect_ratio = :equal,
            label = "trajectory",
            xlim = [xmin[1], xmax[1]],
            ylim = [xmin[2], xmax[2]],
            zlim = [xmin[3], xmax[3]],
            camera = (0.075 * i, 20),
        )
    end


    # fain line
    plot!(plt, eachcol(traj.X)..., lw = 1.5, color = :black, label = nothing, alpha = 0.5)

    plt
end

# -------------------------------- simulation -------------------------------- #

function simulation_frame_1dcan(simulation::Simulation, timems, v::Vector; kwargs...)
    @info "ONE DIM CAN SIM PLOT"
    plt = plot(;
        title = "elapsed: $(round(timems)) ms",
        clims = (0.0, 0.6),
        ylim = [0, 5],
        # aspect_ratio=:equal, 
        grid = false,
    )

    heatmap!(simulation.S')
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

    # offsets = map(offset_for_visual, can.offsets)
    # for (i, offset) in enumerate(offsets)
    #     S = simulation.S[:, i]

    #     # get offset position for plotting
    #     offset = offset .* vec(maximum(can.X; dims = 2))

    #     # plot activity heatmap
    #     x = offset[1] .+ x̄
    #     y = offset[2] .+ ȳ
    #     contourf!(x, y, reshape(S, can.n)', levels = 3)
    # end

    # plot the sum of all activations
    contourf!(x̄, ȳ, reshape(s̄, can.n)', levels = 3)

    # plot input vector direction 
    # x0, y0 = maximum(can.X; dims = 2) ./ 2.1
    # plot!([x0, x0 + v̂[1]], [y0, y0 + v̂[2]], lw = 6, color = :green, alpha=.5, label = nothing)
    # scatter!([x0], [y0], ms = 8, color = :green, label = nothing, alpha=.5, )

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

    th = maximum(s) .* 0.5
    active = s .>= th
    inactive = s .< th

    plt = scatter3d(
        eachrow(M[:, inactive])...,
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
        eachrow(M[:, active])...,
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
    ϕ(x) = moving_average(x, 31)
    if simulation.can.name == "sphere"
        traj = Plots.plot(tj, framen; xmin = [-1.5, -1.5, -1.5], xmax = [1.5, 1.5, 1.5])
        scatter3d!(traj, [x[1]], [x[2]], [x[3]], ms = 5, color = :black, label = "actual")

        # plot decoded trajectory
        plot3d!(
            traj,
            X̄[1:framen, 1],
            X̄[1:framen, 2],
            X̄[1:framen, 3],
            color = :red,
            label = nothing,
            alpha = 1,
            lw = 2,
        )
        scatter3d!(
            traj,
            [X̄[framen, 1]],
            [X̄[framen, 2]],
            [X̄[framen, 3]],
            ms = 7,
            color = :red,
            label = "decoded",
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
    else
        traj = Plots.plot(tj, framen)

        # plot decoded trajectory
        framen > (100 + 2) && plot!(
            traj,
            ϕ(X̄[100:framen, 1]),
            ϕ(X̄[100:framen, 2]),
            color = :red,
            label = nothing,
            alpha = 0.6,
        )
        scatter!(
            traj,
            [X̄[framen, 1]],
            [X̄[framen, 2]],
            ms = 7,
            color = :red,
            label = "decoded",
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
        aspect_ratio = :equal,
        title = "Decoded trajectory",
    )

    scatter!([[x] for x in trajectory.X[1, :]]..., ms = 5, color = :black, label = nothing)

    if d == 2
        plot!(eachcol(X̄)..., lw = 3, color = :red, label = "decoded")
    elseif d == 3
        plot3d!(eachcol(X̄)..., lw = 3, color = :red, label = "decoded")
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
