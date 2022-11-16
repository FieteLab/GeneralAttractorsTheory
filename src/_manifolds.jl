using ForwardDiff: jacobian
using LinearAlgebra



# ---------------------------------------------------------------------------- #
#                                    MOBIUS                                    #
# ---------------------------------------------------------------------------- #
# -------------------- vector fields on the mobius domain -------------------- #
ψ_t(t, θ) = [0, 1]
ψ_t(p) = ψ_t(p...)

ψ_θ1(t, θ) = [(1-sin(4π*t+π/2)) /2 * sign(t), 0] 
ψ_θ1(p) = ψ_θ1(p...)

ψ_θ2(t, θ) = [0, 0]
ψ_θ2(p) = ψ_θ2(p...)



# ---------------------------------------------------------------------------- #
#                                    SPHERE                                    #
# ---------------------------------------------------------------------------- #


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

""" rotation around X axis """
ψx(x, y, z) = (z * ∂y - y * ∂z)
ψx(p) = ψx(p...)

""" rotation around Y axis """
ψy(x, y, z) = (z * ∂x - x * ∂z)
ψy(p) = ψy(p...)

""" rotation around Z axis """
ψz(x, y, z) = (x * ∂y - y * ∂x)
ψz(p) = ψz(p...)

φ = sphere_embedding

""" vector in sphere domain corresponding to vec given by ψx in ℝ³ """
ψxS²(p) = jacobian(φ, p)' * ψx(φ(p))
ψyS²(p) = jacobian(φ, p)' * ψy(φ(p))
ψzS²(p) = jacobian(φ, p)' * ψz(φ(p))



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
