# ---------------------------------------------------------------------------- #
#                                      viz                                     #
# ---------------------------------------------------------------------------- #

function simulation_frame_1dcan(simulation::Simulation, timems, v::Vector; kwargs...)
    plt = plot(;
        title = "elapsed: $(round(timems)) ms",
        grid = false,
    )

    plot!(simulation.g, label="g")
    plot!(simulation.H[:, 1], label="h - pos")
    plot!(simulation.H[:, 2], label="h - neg")
    plt
end


function simulation_frame_2dcan(simulation::Simulation, timems, v::Vector; vmin=0.0, kwargs...)
        main_plot = plot(;
        title = "elapsed: $(round(timems)) ms",
        clims = (vmin, 1.0),
        colorbar=nothing,
        aspect_ratio = :equal,
        grid = false,
    )

    G = simulation.can.G
    g = reshape(simulation.g, G.n)'
    contourf!(g, levels = 3)
    
    m = 1:size(simulation.H, 1)
    plot!(vec(sum(g; dims=1)), lw=2, color="red", label="x")
    plot!(vec(sum(g; dims=2)),lw=2, color="green", label="y")

    m = 1:size(simulation.H, 1)
    ymax = maximum(simulation.H) + 1
    ymin = minimum(simulation.H) - 0.5

    x_shifts = plot(grid=false, title="H₁", ylim=[ymin, ymax])
    plot!(m, simulation.H[:, 1], lw=2, color=:red, label="x-pos")
    plot!(m, simulation.H[:, 2], lw=2, color=:black, label="x-neg")

    y_shifts = plot(grid=false, title="H₂", xlim=[ymin, ymax])
    m = 1:size(simulation.H, 1)
    plot!(simulation.H[:, 3], m, lw=2, color=:red, label="y-pos")
    plot!(simulation.H[:, 4], m, lw=2, color=:black, label="y-neg")

    polar_plot = plot(
        [0, v[1]], [0, v[2]], lw=5, color=:black, label=nothing,
        xlim=(-1.1, 1.1), ylim=(-1.1, 1.1), aspect_ratio=:equal
    )
    scatter!([0],  [0], ms=8, color=:black, label=nothing)

    plot(polar_plot, x_shifts, y_shifts, main_plot; 
        layout = grid(2, 2, heights=[0.3, 0.7], widths=[0.3, 0.7]),
        size=(800, 800),
    )
end



Plots.plot(simulation::Simulation, timems, v::Vector; kwargs...) =
    if simulation.can.d == 1
        simulation_frame_1dcan(simulation, timems, v; kwargs...)
    elseif simulation.can.d == 2
        simulation_frame_2dcan(simulation, timems, v; kwargs...)
    else
        error("Simulation plot for d>2 not implemented")
    end


# ---------------------------------------------------------------------------- #
#                                     utils                                    #
# ---------------------------------------------------------------------------- #
reset(simulation::Simulation) = begin
    simulation.Ṡ .*= 0.0
end
reset(Ṡ::Vector) = Ṡ * 0.0
