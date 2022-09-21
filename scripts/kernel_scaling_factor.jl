import Base.Iterators: product  # cartesian product
import ForwardDiff: jacobian
import LinearAlgebra: norm
import LinearAlgebra: cross as ×
using GLMakie


GLMakie.activate!()
GLMakie.inline!(false)



mfld = :sphere


# ----------------------------------- torus ---------------------------------- #

# DEFINE manifold embedding function φ
if mfld  == :torus
    """ embedding map M → N ⊆ ℝ³  for the Torus"""
    function φ(p::Vector)::Vector
        R, r = 1.0, 0.25
        x, y = p
        return [
            (R + r*cos(x))*cos(y),
            (R + r*cos(x))*sin(y),
            r*sin(x)
        ]
    end
    x =  range(-π, π, length=100)
    y =  range(-π, π, length=50)
    xname, yname = "ϕ₁", "ϕ₂"
else
    """ embedding of the sphere """
    function φ(p::Vector)::Vector
        lon, lat = p

        ls = atan(tan(lat))    # lambda

        return [ cos(ls) * cos(lon),
                2 * cos(ls) * sin(lon),
                sin(ls),
        ]

    end
    x =  range(-π, π, length=100)
    y =  range(-π/2, π/2, length=50)
    xname, yname = "long", "lat"
end


""" partial derivatives of the embedding function """
∂φ∂x(p) = jacobian(φ, p)[:, 1]
∂φ∂y(p) = jacobian(φ, p)[:, 2]


# ---------------------------------- compute --------------------------------- #
# get manifold points
M = hcat([[p...] for p in collect(product(x,y))]...)

# compute product of partials
α(p) = ∂φ∂x(p) × ∂φ∂y(p) |> norm
W = [
    α([x, y]) for x in x, y in y
]

# coordinates of embedded manifold
pts = [
    φ([x, y]) for x in x, y in y
]
X = [p[1] for p in pts]
Y = [p[2] for p in pts]
Z = [p[3] for p in pts]

# --------------------------------- visualize -------------------------------- #
fig = Figure(resolution=(1200, 1200), fontsize=22)
ax = LScene(fig[1, 1], show_axis=true)
pltobj = surface!(
    ax,
    X, Y, Z;
    shading=true,
    ambient=Vec3f(0.65, 0.65, 0.65),
    backlight=1.0f0,
    color=W,
    colormap=:inferno,
    transparency=false,
)
wireframe!(ax, X, Y, Z; transparency=false, color=:black, linewidth=0.5)
cam = cameracontrols(ax.scene)
cam.lookat[] = [0, 0, 200] ./ 1000
cam.eyeposition[] = [5000, 2000, 2000] ./ 1000
cam.upvector[] = [0, 1, 0]
update_cam!(ax.scene, cam)
zoom!(ax.scene, cameracontrols(ax.scene), 1.1)
display(fig)

