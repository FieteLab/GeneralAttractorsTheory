using Plots
import Base.Iterators: product  # cartesian product
import ForwardDiff: jacobian
import LinearAlgebra: norm
import LinearAlgebra: cross as ×

mfld = :torus


# ----------------------------------- torus ---------------------------------- #

# DEFINE manifold embedding function φ
if mfld  == :torus
    """ embedding map M → N ⊆ ℝ³  for the Torus"""
    function φ(p::Vector)::Vector
        R, r = 1.0, 0.75
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
                cos(ls) * sin(lon),
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
W = map(
    p -> ∂φ∂x(p) × ∂φ∂y(p) |> norm, eachcol(M)
)

# --------------------------------- visualize -------------------------------- #
N = mapslices(φ, M, dims=1)

p1 = Plots.contourf(x, y, W, clims=(0, 1.5), 
    title = "$mfld: ||∂φ∂xₚ × ∂φ∂yₚ||",
    xlabel=xname, ylabel=yname
) 


#? 3D scatter plot
p2 = Plots.scatter(
    N[1, :], N[2, :], N[3, :], marker_z=W, msw=0, ms=10, 
    camera=(0, 90),
    xlim=[-1.2, 1.2], ylim=[-1.2, 1.2], zlim=[-1.2, 1.2],
    title = "Embedded $mfld top view"
)
display(Plots.plot(p1, p2, size=(800, 600), layout=grid(2,1)))

# TODO fix plotting bug
# TODO think of what S meant