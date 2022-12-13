"""
Collection of code useful to visualize manifolds
"""
module ManifoldUtils
using Distances

import ..GeneralAttractors: SphericalDistance, MobiusEuclidean, lerp
import ..GeneralAttractors: sphere_embedding, mobius_embedding

export AbstractManifold, CoverSpace
export ℝ², T, S²

include("_manifolds.jl")


# ---------------------------------------------------------------------------- #
#                                 VECTOR FIELDS                                #
# ---------------------------------------------------------------------------- #

abstract type AbstractVectorField end

struct ConstantVectorField <: AbstractVectorField
    d::Int
    i::Int
end

function (ψ::ConstantVectorField)(::Vector)
    v = zeros(ψ.d)
    v[ψ.i] = 1
    return v
end

struct VectorField <: AbstractVectorField
    f::Function
end
(ψ::VectorField)(x::Vector) = ψ.f(x)



# ---------------------------------------------------------------------------- #
#                                   MANIFOLDS                                  #
# ---------------------------------------------------------------------------- #
abstract type AbstractManifold end


function Base.rand(m::AbstractManifold)
    d = length(m.xmin)
    x = zeros(d)
    for i = 1:d
        x[1] = rand(m.xmin[i]:0.025:m.xmax[i])
    end
    x
end

""" Ensure x in M """
apply_boundary_conditions!(x, ::AbstractManifold) = (x, ones(length(x)))

# ---------------------------------------------------------------------------- #
#                                     RING                                     #
# ---------------------------------------------------------------------------- #
struct Ring <: AbstractManifold
    xmin::Vector
    xmax::Vector
    ψs::Vector{AbstractVectorField}
    metric::Metric
end
Ring() = Ring([0], [2π], 
    [VectorField(ring_ψ)], 
    # [ConstantVectorField(1, 1)],
    PeriodicEuclidean([2π]))

apply_boundary_conditions!(x, ::Ring) = mod.(x, 2π), 1


# ---------------------------------------------------------------------------- #
#                                     PLANE                                    #
# ---------------------------------------------------------------------------- #
struct Manifoldℝ² <: AbstractManifold
    xmin::Vector
    xmax::Vector
    ψs::Vector{AbstractVectorField}
    metric::Metric
end
Manifoldℝ²(m) = Manifoldℝ²(
    [-m, -m],
    [m, m],
    # [ConstantVectorField(2, 1), ConstantVectorField(2, 2)],
    [VectorField(ℝ²_ψ1), VectorField(ℝ²_ψ2)],
    Euclidean(),
)
ℝ² = Manifoldℝ²(100)


# ---------------------------------------------------------------------------- #
#                                     TORUS                                    #
# ---------------------------------------------------------------------------- #
struct Torus <: AbstractManifold
    xmin::Vector
    xmax::Vector
    ψs::Vector{AbstractVectorField}
    metric::Metric
end
Torus() = Torus(
    [0, 2π],
    [0, 2π],
    [ConstantVectorField(2, 1), ConstantVectorField(2, 2)],
    PeriodicEuclidean([2π, 2π]),
)
T = Torus()

function apply_boundary_conditions!(x::Vector, ::Torus)
    return mod.(x, 2π), ones(length(x))
end


# ---------------------------------------------------------------------------- #
#                                    SPHERE                                    #
# ---------------------------------------------------------------------------- #
struct Sphere <: AbstractManifold
    xmin::Vector
    xmax::Vector
    ψs::Vector{AbstractVectorField}
    metric::Metric
end
S² = Sphere(
    [-π, -π / 2],
    [π, π / 2],
    [VectorField(ψx), VectorField(ψy), VectorField(ψz)],
    SphericalDistance(),
)

# ---------------------------------------------------------------------------- #
#                                    MOBIUS                                    #
# ---------------------------------------------------------------------------- #
struct Mobius <: AbstractManifold
    xmin::Vector
    xmax::Vector
    ψs::Vector{AbstractVectorField}
    metric::Metric
end
Mobius() = Mobius(
    [-0.75, 0],
    [ 0.75, 2π],
    [VectorField(ψ_t), VectorField(ψ_θ1), VectorField(ψ_θ2)],
    MobiusEuclidean(),
)


""" 
Correct the position vector x to ensure
that it's on the MB. If it's too much to the side
on the non periodic dimension, put it at the manifold's boundary,
if it's along the periodic dimension gets its position module 2π.
"""
function apply_boundary_conditions!(x::Vector, m::Mobius)
    vel_correction_facors = [1, 1]
    # non periodic dimension
    δ = 0.2  # padding around boundary to account for bump size
    if x[1] <= m.xmin[1] + δ
        x[1] = m.xmin[1] + δ
        vel_correction_facors[1] = 0
    elseif x[1] >= m.xmax[1] - δ
        x[1] = m.xmax[1] - δ
        vel_correction_facors[1] = 0
    end

    # periodic dimension
    if x[2] >= 2π
        x[2] = x[2] - 2π
        x[1] = -x[1]
    elseif x[2] <= 0
        x[2] = 2π + x[2]
        x[1] = -x[1]
    end
    return x, vel_correction_facors
end


# ---------------------------------------------------------------------------- #
#                                 COVER SPACES                                 #
# ---------------------------------------------------------------------------- #

struct CoverSpace
    M::AbstractManifold  # variable manifold (cover space)
    N::AbstractManifold  # neural manifold (covered)
    ρ::Function          # cover map | it's ρ (\rho) not 'p'
    ρⁱ::Function         # inverse of the cover map
end

CoverSpace(M, N) = CoverSpace(M, N, identity, identity)
CoverSpace(M) = CoverSpace(M, M)

identity(x) = x


end
