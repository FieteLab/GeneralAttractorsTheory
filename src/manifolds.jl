"""
Collection of code useful to visualize manifolds
"""
module Manifolds

    export AbstractManifold
    export ℝ²

    abstract type AbstractManifold end


    struct Manifoldℝ² <: AbstractManifold
    end

    Base.string(::Manifoldℝ²) = "ℝ²"

    ℝ² = Manifoldℝ²()

end
