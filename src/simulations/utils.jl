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
        title = "M: $(traj.M)",
        label = "trajectory",
        grid = false,
        aspect_ratio = :equal,
    )
end


function Plots.plot(traj::Trajectory, i::Int)
    d = size(traj.X, 2)
    d != 2 && error("not implemented")

    plot(
        traj.X[1:i, 1],
        traj.X[1:i, 2],
        lw = 3,
        color = :black,
        title = "M: $(traj.M)",
        label = "trajectory",
        grid = false,
        aspect_ratio = :equal,
    )
end

function simulation_frame_1dcan(simulation::Simulation, timems, v::Vector; kwargs...)
    @info "ONE DIM CAN SIM PLOT"
    can = simulation.can
    plt = plot(;
        title = "elapsed: $(round(timems)) ms",
        clims = (0.0, 0.6),
        ylim = [0, 5],
        # aspect_ratio=:equal, 
        grid = false,
    )

    heatmap!(simulation.S')
    # hline!([1.5], lw=5, color=:white, label=nothing)

    # x0 = size(simulation.S, 1)/2
    # v̂ = v ./ v[1]
    # plot!(
    #     [x0, x0+v̂[1]*x0/4], [3, 3], lw=10, color=:black, label=nothing
    # )
    # scatter!([x0], [3], ms=10, color=:black, label=nothing)

    plt
end


function simulation_frame_2dcan(simulation::Simulation, timems, v::Vector; kwargs...)
    can = simulation.can
    plt = plot(;
        title = "elapsed: $(round(timems)) ms",
        clims = (0.0, 1.0),
        aspect_ratio = :equal,
        grid = false,
        size = simulation.can.n .* 10,
    )

    h = maximum(simulation.can.X) / 2.5
    v̂ = v ./ norm(v) .* h
    for i = 1:can.d*2
        S = simulation.S[:, i]

        # get offset position for plotting
        offset = can.A[i, :] ./ maximum(can.A) .* vec(maximum(simulation.can.X; dims = 2))

        # plot activity heatmap
        x = offset[1] .+ range(0, maximum(simulation.can.X[1, :]), length = can.n[1])
        y = offset[2] .+ range(0, maximum(simulation.can.X[2, :]), length = can.n[2])
        contourf!(x, y, reshape(S, can.n)', levels = 3)

        # plot input vector direction 
        x0, y0 = maximum(simulation.can.X; dims = 2) ./ 2.1

        plot!([x0, x0 + v̂[1]], [y0, y0 + v̂[2]], lw = 6, color = :green, label = nothing)
        scatter!([x0], [y0], ms = 8, color = :green, label = nothing)

        # # add text with Aᵢ
        # Ai = string(round.(ai; digits=3))
        # vi = ai ⋅ v

        # annotate!(
        #     [offset[1]+5, x0], 
        #     [offset[2]+20, y0+20],      
        #     [
        #         text(" I: $i\nAi: $Ai\nAv: $(round(vi; digits=4))", :white, :left, 8),
        #         text("v⃗: $(round.(v; digits=2))", :black, :center, 8),
        #     ]
        # )
    end
    plt
end



function Plots.plot(simulation::Simulation, timems, framen, v::Vector; kwargs...)
    if simulation.can.d == 1
        pop_activity = simulation_frame_1dcan(simulation, timems, v; kwargs...)
    elseif simulation.can.d == 2
        pop_activity = simulation_frame_2dcan(simulation, timems, v; kwargs...)
    else
        error("Simulation plot for d>2 not implemented")
    end

    traj = Plots.plot(traj, framen)


    plot(traj, pop_activity, size = (1000, 800))
end


# ---------------------------------------------------------------------------- #
#                                     utils                                    #
# ---------------------------------------------------------------------------- #
reset(simulation::Simulation) = begin
    simulation.Ṡ .*= 0.0
end
reset(Ṡ::Vector) = Ṡ * 0.0
