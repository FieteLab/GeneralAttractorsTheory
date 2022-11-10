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
ℝ² = Manifoldℝ²([-100, 100], [100, 100])

struct Torus <: AbstractManifold
    xmin::Vector
    xmax::Vector
end
T = Torus([-32, -32], [32, 32])


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
    ρ::Function          # cover map
end


end
