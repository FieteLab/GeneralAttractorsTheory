using ForwardDiff: jacobian
using LinearAlgebra
using Rotations


R = RotXYZ(90, 90, 90)


"""
Classic sphere embedding in ℝ³
"""
function sphere_embedding end

sphere_embedding(p) = sphere_embedding(p...)
function sphere_embedding(lon, lat) 
    ls = atan(tan(lat)) 
    return [ 
            cos(ls) * cos(lon),
            cos(ls) * sin(lon),
            sin(ls),
    ]
end


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
function metric_deformation(φ::Function, p::AbstractVector)::Tuple{Number, Number}
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