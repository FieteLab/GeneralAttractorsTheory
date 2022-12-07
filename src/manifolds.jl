"""
Collection of code useful to visualize manifolds
"""
module ManifoldUtils
using Distances

import ..GeneralAttractors: SphericalDistance, MobiusEuclidean
import ..GeneralAttractors: sphere_embedding, mobius_embedding

export AbstractManifold, CoverSpace
export ℝ², T, S²

include("_manifolds.jl")


# ---------------------------------------------------------------------------- #
#                                   MANIFOLDS                                  #
# ---------------------------------------------------------------------------- #
abstract type AbstractManifold end

struct Ring <: AbstractManifold
    xmin::Vector
    xmax::Vector
    metric::Metric
end
Ring() = Ring([0,], [2π, ], PeriodicEuclidean([2π]))

struct Manifoldℝ² <: AbstractManifold
    xmin::Vector
    xmax::Vector
    metric::Metric
end
Manifoldℝ²(m) = Manifoldℝ²([-m, -m], [m, m], Euclidean())
ℝ² = Manifoldℝ²(250)

struct Torus <: AbstractManifold
    xmin::Vector
    xmax::Vector
    metric::Metric
end
Torus() = Torus([0, 2π], [0,  2π], PeriodicEuclidean([2π, 2π]))
T = Torus()


struct Sphere <: AbstractManifold
    xmin::Vector
    xmax::Vector
    metric::Metric
end
S² = Sphere([-π, -π / 2], [π, π / 2], SphericalDistance())


struct Mobius <: AbstractManifold
    xmin::Vector
    xmax::Vector
    metric::Metric
end
Mobius() = Mobius([-1 / 2, 1 / 2], [0, 2π], MobiusEuclidean())

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
