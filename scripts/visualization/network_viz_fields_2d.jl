using GeneralAttractors.Simulations
import MyterialColors: indigo, salmon_dark, black

using GeneralAttractors.ManifoldUtils: Mobius, ψ_t, ψ_θ1, ψ_θ2
using Plots

"""
Visualization of the vector fields on the mobiuscan
"""

include("../networks/torus.jl")

can = toruscan
scaling = 1
colors = [black, black, indigo, indigo, salmon_dark, salmon_dark]



p = plot(
    aspect_ratio = :equal,
    size = (600, 600),
    grid = false,
)


n = size(can.X, 2)
for i in 1:4:n
    x = can.X[:, i]
    scatter!([[x] for x in x]..., label=nothing, color=:black, ms=3)

    for (j, o) in enumerate(can.Ω)
        j % 2 == 0 && continue
        v = o(x)
        plot!(
            [x[1], x[1]+v[1]*scaling],
            [x[2], x[2]+v[2]*scaling],
            lw=2, color=colors[j], label=nothing
        )
    end
end

# for t = -1/2:0.4:1/2, θ = 0:0.25:2π
#     d1 = ψ1([t, θ]) .* scaling
#     plot!([t, t + d1[1]], [θ, θ + d1[2]], color = :black, label = nothing)
#     plot!([t + 1.5, t + d1[1] + 1.5], [θ, θ + d1[2]], color = :black, label = nothing)

#     d1 = ψ2([t, θ]) .* scaling
#     plot!([t, t + d1[1]], [θ, θ + d1[2]], color = :red, label = nothing)

#     d1 = ψ3([t, θ]) .* scaling
#     plot!([t + 1.5, t + d1[1] + 1.5], [θ, θ + d1[2]], color = :green, label = nothing)
# end

p
