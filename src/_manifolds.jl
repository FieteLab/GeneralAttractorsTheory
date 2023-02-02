using ForwardDiff: jacobian
using LinearAlgebra
using Distances: euclidean

# ---------------------------------------------------------------------------- #
#                                     RING                                     #
# ---------------------------------------------------------------------------- #
ring_ψ(x) = ring_ψ(x...)
ring_ψ(x::Number) = (sin(x)+1.1)/2

# ---------------------------------------------------------------------------- #
#                                  PLANE/TORUS                                 #
# ---------------------------------------------------------------------------- #

# curl free vector field
torus_ψ1(x, y) = [cos(x)/2 + 1, 0]
torus_ψ2(x, y) = [0, cos(y)/2 + 1 ]

# curly vector field
# torus_ψ1(x, y) = [
#     cos(x),
#     sin(x),
#     ]
# torus_ψ2(x, y) =  [[0, 1] [-1, 0]]  * torus_ψ1(x, y)


torus_ψ1(x) = torus_ψ1(x...)
torus_ψ2(x) = torus_ψ2(x...)


# ---------------------------------------------------------------------------- #
#                                    MOBIUS                                    #
# ---------------------------------------------------------------------------- #

""" constant """
MB_ψ1(t, θ) = [0, 1]
MB_ψ2(t, θ) = [1, 0]
MB_ψ3(t, θ) = [-1, 0]

MB_ψ1(p) = MB_ψ1(p...)
MB_ψ2(p) = MB_ψ2(p...)
MB_ψ3(p) = MB_ψ3(p...)


# ---------------------------------------------------------------------------- #
#                                    SPHERE                                    #
# ---------------------------------------------------------------------------- #

"""
(almost) equally spaced points on the sphere.
https://stackoverflow.com/questions/9600801/evenly-distributing-n-points-on-a-sphere
"""
function fibonacci_sphere(n = 1000)
    points = zeros(3, n)
    ϕ = π * (3 - √5)  # golden angle in radians

    for i = 1:n
        y = 1 - (i / float(n - 1)) * 2  # y goes from 1 to -1
        radius = √(Complex(1 - y * y)) |> real  # radius at y

        θ = ϕ * i  # golden angle increment

        x = cos(θ) * radius
        z = sin(θ) * radius

        points[:, i] = [x, y, z]
    end
    return points
end

"""
Killing fields on the unit sphere

These are defined as rotational vector fields
on the unit sphere in ℝ³ with rotations around X, Y, Z.
These are killing fields and are divergence free.


Given a point p=(x, y, z) ∈ S² ⊂ ℝ²,
get a vector (fx(p)∂x, fy(p)∂y, fz(p)∂z)
tangent to the sphere and correpsonding to a rotation.
"""

∂x = [1, 0, 0]
∂y = [0, 1, 0]
∂z = [0, 0, 1]


normalize(x) = norm(x) > 0 ? x ./ norm(x) : x


""" rotation around X axis """
# sphere_ψx(x, y, z)::Vector = (z * ∂y - y * ∂z) |> normalize
sphere_ψx(x, y, z)::Vector = [0, -z, y] # |> normalize
sphere_ψx(p) = sphere_ψx(p...)

""" rotation around Y axis """
# sphere_ψy(x, y, z)::Vector = (z * ∂x - x * ∂z) |> normalize
sphere_ψy(x, y, z)::Vector = [z, 0, -x] # |> normalize
sphere_ψy(p) = sphere_ψy(p...)

""" rotation around Z axis """
# sphere_ψz(x, y, z)::Vector = (x * ∂y - y * ∂x) |> normalize
sphere_ψz(x, y, z)::Vector = [-y, x, 0] # |> normalize
sphere_ψz(p) = sphere_ψz(p...)



# ---------------------------------------------------------------------------- #
#                                   DIFF GEOM                                  #
# ---------------------------------------------------------------------------- #
""" 
First fundamental form of an embedding φ: M → N at a point p ∈ M
"""
function 𝐈(φ::Function, p::AbstractVector)::Matrix
    J = jacobian(φ, p)
    return J' * J
end

"""
    metric_deformation(φ::Function, p::AbstractVector)::Tuple{Number, Number}

Eigenvalues of the first fundamental form
"""
function metric_deformation(φ::Function, p::AbstractVector)::Tuple{Number,Number}
    λ₁, λ₂ = eigen(𝐈(φ, p)).values
    return λ₁, λ₂
end

"""
    area_deformation(φ::Function, p::AbstractVector)::Float64

Area deformation given an embedding φ.
"""
function area_deformation(φ::Function, p::AbstractVector)::Float64
    λ₁, λ₂ = metric_deformation(φ, p)
    λ₁ * λ₂
end
