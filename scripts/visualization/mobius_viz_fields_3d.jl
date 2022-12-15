using GeneralAttractors.Simulations

using GeneralAttractors.ManifoldUtils: Mobius, mobius_embedding
using Plots

"""
Visualization of the vector fields on the mobiuscan
"""

include("../networks/mobius.jl")


ψ1, ψ2, ψ3 = mobiuscan.C.M.ψs
scaling = .07



plt = plot(
    aspect_ratio = :equal,
    size = (600, 600),
    grid = false,
    # camera=(90, 0)
    camera=(90, 0)

)


for t = -1/2:0.1:1/2, θ = 0:0.25:(2π-0.25)
    m = [[x] for x in mobius_embedding([t, θ])]
    scatter3d!(m..., color=:black, label=nothing, ms=2)

    for (ψ, c) in zip((ψ1, ψ2, ψ3), (:black, :red, :green))
        p = [t, θ] .+ ψ([t, θ]) .* scaling |> mobius_embedding
        plot3d!(
            [m[1][1], p[1]],
            [m[2][1], p[2]],
            [m[3][1], p[3]],
            lw=4, color=c, label=nothing
        )
    end
end

plt