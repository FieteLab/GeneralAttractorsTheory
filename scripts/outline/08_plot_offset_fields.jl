include("settings.jl")
using MyterialColors
import Base.Iterators: product as ×


for (name, can) in cans
    name == "sphere" && begin
        @warn "Skipping sphere"
        continue
    end
    xmin, xmax = can.C.N.xmin, can.C.N.xmax

    # get 100 points between xmin and xmax
    Δ = (xmax .- xmin)/5
    pts = (range(xmin[1], xmax[1], step=Δ[1]) × range(xmin[2], xmax[2], step=Δ[2])) |> collect
    second(x) = x[2]
    pts = hcat(vcat(first.(pts)...), vcat(second.(pts)...))

    # plts = []
    plt = plot(axis=:off, size=(600,600), title=name, aspect_ratio=:equal, grid=false)
    plot!([(xmin[1], xmin[2]), (xmax[1], xmin[2]), (xmax[1], xmax[2]), (xmin[1], xmax[2])],
    seriestype=:shape,
    fillcolor="#E0DCE8",
    linecolor=:black,
    linewidth=0.6, label=nothing
    )

    colors = [salmon, salmon_darker, indigo, indigo_darker, teal, teal_darker]
    for (j, (o, c)) in enumerate(zip(can.offsets, colors))
        j % 2 != 1 && continue
        for pt in eachrow(pts)
            # draw "offset"
            v = pt + (o(pt)/5)
            plot!(
                [pt[1], v[1]], [pt[2], v[2]],
                linecolor=c,
                linewidth=3, label=nothing

            )
        end
    end

    scatter!(pts[:, 1], pts[:, 2], markersize=4.5, color="black", label=nothing)
    display(plt)
end