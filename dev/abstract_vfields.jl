using Plots
using GeneralAttractors
import GeneralAttractors.ManifoldUtils: sphere_embedding, ψx, ψy, ψz, fibonacci_sphere
using MyterialColors
import GeneralAttractors: by_column


x = -π:.1:π
y = -π/:.1:π/2


# M = hcat([[x, y] for x in x for y in y])
# N = hcat(map(r -> sphere_embedding(r), eachrow(M))...)
N = fibonacci_sphere(500)


plots = []
for (vfield, color) in zip((ψx, ψy, ψz),(blue_grey_dark, salmon_dark, indigo_dark))
    plt = plot(grid=false,
    xticks=nothing, yticks=nothing, zticks=nothing,
    camera=(90, 00),
    xlim=[-1, 1],
    ylim=[-1, 1],
    zlim=[-1, 1],
)

s = scatter3d
    for p in eachcol(N[:, 1:2:end])
        p[1] < 0 && continue
        vx = vfield(p) .* .2

        plot3d!(
            [p[1], p[1]+vx[1]],
            [p[2], p[2]+vx[2]],
            [p[3], p[3]+vx[3]],
            lw=2, color=color, label=nothing
        )
        scatter3d!(
            [p[1],],
            [p[2],],
            [p[3],],
            ms=2,
            lw=2, color=color, label=nothing,
            msa=0, msw=0
        )

    end
    push!(plots, plt)
end

s =  scatter3d(eachrow(N[:, N[1, :] .> 0])..., label=nothing, grid=false,
    xticks=nothing, yticks=nothing, zticks=nothing,
    camera=(90, 00),
    xlim=[-1, 1],
    ylim=[-1, 1],
    zlim=[-1, 1],
    ms=2,
    lw=2, color=:black,
    msa=0, msw=0
)


plot(plots...,s, layout=(1, 4), size=(800, 400))


