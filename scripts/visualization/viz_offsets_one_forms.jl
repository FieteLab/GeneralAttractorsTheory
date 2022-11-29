using Plots

using GeneralAttractors


# include("../networks/torus.jl")


can = toruscan

i = 1
scale = .5
δx = 1
X = can.C.N.xmin[1]:δx:can.C.N.xmax[2] |> collect
Y = can.C.N.xmin[1]:δx:can.C.N.xmax[2]
pts = hcat([[x, y] for x in X for y in Y]...)

offset_plot = plot(grid=false, title="offsets")
oneform_plot = plot(grid=false, title="one form")

scatter!(offset_plot, eachrow(pts)..., label=nothing, color=:black, ms=2)
scatter!(oneform_plot, eachrow(pts)..., label=nothing, color=:black, ms=2)

for (x, y) in eachcol(pts)
    o = can.offsets[i](x, y) .* scale
    plot!(
        offset_plot, [x, x+o[1]],[y, y+o[2]], lw=2, color=:green, label=nothing
    )

    ω = can.Ω[i]([x, y]) .* scale
    plot!(
        oneform_plot, [x, x+ω[1]],[y, y+ω[2]], lw=2, color=:red, label=nothing
    )
end


plot(offset_plot, oneform_plot)