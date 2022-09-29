import Base.Iterators: product  # cartesian product
import ForwardDiff: jacobian
import LinearAlgebra: norm, eigen
import LinearAlgebra: cross as ×
using GLMakie
using GeneralAttractors
using GeneralAttractors.ManifoldUtils


GLMakie.inline!(true)

mfld = :torus


# ----------------------------------- torus ---------------------------------- #

# DEFINE manifold embedding function φ
if mfld  == :torus
    embedding = TorusEmbedding()
    x =  range(-π, π, length=25)
    y =  range(-π, π, length=40)
    xname, yname = "ϕ₁", "ϕ₂"
elseif mfld == :sphere
    embedding = SphereEmbedding()
    x =  range(-π, π, length=100)
    y =  range(-π/2, π/2, length=50)
    xname, yname = "long", "lat"

elseif mfld == :mobius
    embedding = MobiusEmbedding()
    x =  range(-π, π, length=50)
    y =  range(-1, 1, length=10)
    xname, yname = "x", "y"

elseif mfld == :cylinder
    embedding = CylinderEmbedding()
    x =  range(-π, π, length=100)
    y =  range(-1, 1, length=5)
    xname, yname = "x", "y"

else
    throw("mfld not recognized $mfld")
end


""" partial derivatives of the embedding function """
∂φ∂x(p) = jacobian(embedding.φ, p)[:, 1]
∂φ∂y(p) = jacobian(embedding.φ, p)[:, 2]



# ---------------------------------- compute --------------------------------- #
# get manifold points
M = hcat([[p...] for p in collect(product(x,y))]...)

# equation for the normal of the surface at a point p
n(p) = begin
    a = ∂φ∂x(p) × ∂φ∂y(p)
    a ./ norm(a)
end 
    
# compute product of partials
# α(p) = ∂φ∂x(p) × ∂φ∂y(p) |> norm

# compute the eigenvalues of the first intrinsic form
function α(p)
    J = jacobian(embedding.φ, p)
    I = J' * J
    λ₁, λ₂ = eigen(I).values
    λ₁*λ₂
end

W = [
    α([x, y]) for x in x, y in y
]


visualize_manifold(
    x, y, embedding; color=W, cmap=:inferno, transparency=false
)