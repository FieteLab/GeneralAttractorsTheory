using GLMakie

using GeneralAttractors
using GeneralAttractors.Simulations
using GeneralAttractors.ManifoldUtils
import GeneralAttractors.ManifoldUtils: sphere_embedding

"""
Animate a simulation given it's history `h`

TODO: load simulation and h so that we don't have to run first
TODO: save in correct place/name
"""

# get embedded manifold coordinates
coords = map(sphere_embedding, eachcol(simulation.can.X))

X = first.(coords)
Y = [c[2] for c in coords]
Z = last.(coords)

# get activity
m, _, n = size(h.S)
S = reshape(sum(h.S, dims = 2), (m, n))

# plot
s = Observable(S[:, 1])
fig = Figure(resolution = (1200, 1200), viewmode = :fitzoom)
ax = LScene(
    fig[1, 1],
    scenewk = (; padding = (0, 0, 0)),
    # scenekw = (; limits=Rect3f(Vec3f(-1, -1, -1),Vec3f(2, 2, 2)))
)

pltobj = GLMakie.scatter!(
    ax,
    X,
    Y,
    Z;
    shading = false,
    markersize = 150,
    strokewidth = 0,
    color = s,
    fxaa = false,
    ssao = false,
)

record(fig, "time_animation.mp4"; framerate = 30) do io
    for t = 1:n
        t % 100 == 0 && println("$t/$n")
        s[] = S[:, t]   # updates pltobj
        recordframe!(io)
    end
end


println("Done")
