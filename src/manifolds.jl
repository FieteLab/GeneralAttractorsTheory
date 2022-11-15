"""
Collection of code useful to visualize manifolds
"""
module ManifoldUtils

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
end
Manifoldℝ²(m) = Manifoldℝ²([-m, -m], [m, m])
ℝ² = Manifoldℝ²([-250, -250], [250, 250])

struct Torus <: AbstractManifold
    xmin::Vector
    xmax::Vector
end
Torus(m) = Torus([-m/2, -m], [m/2, m/2])
T = Torus(32)


struct Sphere <: AbstractManifold
    xmin::Vector
    xmax::Vector
end
S² = Sphere([-π, -π / 2], [π, π / 2])

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
