"""
Collection of code useful to visualize manifolds
"""
module ManifoldUtils
using Distances

import ..GeneralAttractors: SphericalDistance, MobiusEuclidean

export AbstractManifold, CoverSpace
export ℝ², T, S²

include("_manifolds.jl")


# ---------------------------------------------------------------------------- #
#                                   MANIFOLDS                                  #
# ---------------------------------------------------------------------------- #
abstract type AbstractManifold end


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
Torus(m) = Torus([-m/2, -m], [m/2, m/2], PeriodicEuclidean([m, m]))
T = Torus(32)


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
Mobius(m) = Mobius([0, 0], [m, m], MobiusEuclidean(m))

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
