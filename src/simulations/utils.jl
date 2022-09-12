# ---------------------------------------------------------------------------- #
#                                      viz                                     #
# ---------------------------------------------------------------------------- #


function Plots.plot(simulation::Simulation, timems, v::Vector; kwargs...)
    can = simulation.can
    plt = plot(;
        title="elapsed: $(round(timems)) ms", 
        # clims=(0.0, 0.05),
        aspect_ratio=:equal, 
        grid=false,
    )

    
    for i in 1:can.d * 2
        S = simulation.S[:, i]
        S = reshape(S, can.n)'

        # get offset position for plotting
        ai = can.offset_directions[i]
        offset = ceil.(ai) .* can.n
        x = collect(1:can.n[2]) .+ offset[1]
        y = collect(1:can.n[1]) .+ offset[2]

        heatmap!(x, y, S)
        x0, y0 = can.n[1]/2, can.n[2]/2
        plot!([x0, x0+v[1]*20], [y0, y0+v[2]*20], lw=6, color=:black, label=nothing)
        scatter!([x0], [y0], ms=8, color=:black, label=nothing)

        # add text with Aᵢ
        Ai = string(round.(ai; digits=3))
        vi = ai ⋅ v
        annotate!(
            [offset[1]+5, x0], 
            [offset[2]+20, y0+20],      
            [
                text(" I: $i\nAi: $Ai\nVi: $(round(vi; digits=4))", :white, :left, 7),
                text("v⃗: $(round.(v; digits=2))", :black, :center, 8),
            ]
        )
    end
    plt
end


# ---------------------------------------------------------------------------- #
#                                     utils                                    #
# ---------------------------------------------------------------------------- #
reset(simulation::Simulation) = begin
    simulation.Ṡ .*= 0.0
end
reset(Ṡ::Vector) = Ṡ * 0.0

