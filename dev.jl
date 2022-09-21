using GeneralAttractors
using GLMakie
using Distances
import Base.Iterators: product  # cartesian product
import ForwardDiff: jacobian
import LinearAlgebra: norm
import LinearAlgebra: cross as ×

GLMakie.activate!()
GLMakie.inline!(true)

function φ(p::Vector)::Vector
    lon, lat = p

    ls = atan(tan(lat))    # lambda

    return [ cos(ls) * cos(lon),
            cos(ls) * sin(lon),
            sin(ls),
    ]

end

xs =  range(-π, π, length=50)
ys =  range(-π/2, π/2, length=50)

# M = hcat([[p...] for p in collect(product(x,y))]...)
# N = mapslices(φ, M, dims=1)

pts = [
    φ([x, y]) for x in xs, y in ys
]
X = [p[1] for p in pts]
Y = [p[2] for p in pts]
Z = [p[3] for p in pts]



fig = Figure(resolution=(1200, 800), fontsize=22)
ax = LScene(fig[1, 1], show_axis=true)
pltobj = surface!(
    ax,
    X, Y, Z;
    shading=true,
    ambient=Vec3f(0.65, 0.65, 0.65),
    backlight=1.0f0,
    # color=sqrt.(X1 .^ 2 .+ Y1 .^ 2 .+ Z1 .^ 2),
    # colormap=Reverse(:bone_1),
    transparency=true,
)
wireframe!(ax, X, Y, Z; transparency=true, color=:gray, linewidth=0.5)
zoom!(ax.scene, cameracontrols(ax.scene), 0.98)
display(fig)


