"""
Collection of code useful to visualize manifolds
"""
module Manifolds

export AbstractManifold, CoverSpace
export ℝ², T




# ---------------------------------------------------------------------------- #
#                                   MANIFOLDS                                  #
# ---------------------------------------------------------------------------- #
abstract type AbstractManifold end


struct Manifoldℝ² <: AbstractManifold end
ℝ² = Manifoldℝ²()

struct Torus <: AbstractManifold end
T = Torus()

# ---------------------------------------------------------------------------- #
#                                 COVER SPACES                                 #
# ---------------------------------------------------------------------------- #

struct CoverSpace
    M::AbstractManifold  # variable manifold (cover space)
    N::AbstractManifold  # neural manifold (covered)
    ρ::Function          # cover map
end


end
