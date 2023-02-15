using Plots

using GeneralAttractors
import GeneralAttractors.ManifoldUtils: fibonacci_sphere

pyplot()

pts = fibonacci_sphere(1000)


plt = scatter3d(eachrow(pts)...)


X(x, y , z) = [0, -z, y]
Y(x, y , z) = [z, 0, -x]
Z(x, y , z) = [-y, x, 0]

function plot_vec!(p, ψ, color)
    v = ψ(p...) .* .25


    plot!(
        [p[1], p[1] + v[1]],
        [p[2], p[2] + v[2]],
        [p[3], p[3] + v[3]],
        color = color,
        lw = 2,
        label = nothing,
        arrow=true,
    )
end


for p in eachcol(pts)
    plot_vec!(p, X, :red)
    # plot_vec!(p, Y, :green)
end


display(plt)