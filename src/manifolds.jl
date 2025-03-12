"""
Collection of code useful to visualize manifolds
"""
module ManifoldUtils
using Distances
import Term.Repr: @with_repr, termshow
using Interpolations
import ..GeneralAttractors: SphericalDistance, MobiusEuclidean, lerp, KleinBottleEuclidean
import ..GeneralAttractors: sphere_embedding, mobius_embedding, klein_embedding

export AbstractManifold, CoverSpace
export ℝ², T, S²
export Ring, Mobius, Sphere, Torus, Manifoldℝ², Cylinder, Line, KleinBottle

include("_manifolds.jl")

export MB_ψ1, MB_ψ2, MB_ψ3
export ring_ψ
export torus_ψ1, torus_ψ2
export sphere_ψx, sphere_ψy, sphere_ψz, fibonacci_sphere
export make_space_filling_trajectory

# ---------------------------------------------------------------------------- #
#                                 VECTOR FIELDS                                #
# ---------------------------------------------------------------------------- #

abstract type AbstractVectorField end

@with_repr struct ConstantVectorField <: AbstractVectorField
    d::Int
    i::Int
end

function (ψ::ConstantVectorField)(::Vector)
    v = zeros(ψ.d)
    v[ψ.i] = 1
    return v
end

@with_repr struct VectorField <: AbstractVectorField
    f::Function
end
(ψ::VectorField)(x::Vector) = ψ.f(x)



# ---------------------------------------------------------------------------- #
#                                   MANIFOLDS                                  #
# ---------------------------------------------------------------------------- #
abstract type AbstractManifold end


function Base.rand(m::AbstractManifold; δ=0)
    d = length(m.xmin)
    x = zeros(d)
    for i = 1:d
        x[i] = rand((m.xmin[i]+δ):0.001:(m.xmax[i]-δ))
    end
    x
end

""" Ensure x in M """
apply_boundary_conditions!(x, ::AbstractManifold) = (x, ones(length(x)))


# ---------------------------------------------------------------------------- #
#                                     LINE                                     #
# ---------------------------------------------------------------------------- #

@with_repr struct Line <: AbstractManifold
    name::String
    xmin::Vector
    xmax::Vector
    ψs::Vector{AbstractVectorField}
    metric::Metric
    d::Int  # dimensionality
    periodic_dimensions::Vector # label for which dimensions are periodic
end

Line(; d=6) = Line(
    "line",
    [-d],
    [d],
    [ConstantVectorField(1, 1)],
    Euclidean(),
    1, [0]
)

function apply_boundary_conditions!(x, m::Line)
    δ = 0.5
    x[1] = max(m.xmin[1]+δ, min(m.xmax[1]-δ, x[1]))
    return x, 1
end

function Base.rand(m::Line; δ=0.25)
    x = rand((m.xmin[1]+δ):0.001:(m.xmax[1]-δ))
    [x,]
end


# ---------------------------------------------------------------------------- #
#                                     RING                                     #
# ---------------------------------------------------------------------------- #
@with_repr struct Ring <: AbstractManifold
    name::String
    xmin::Vector
    xmax::Vector
    ψs::Vector{AbstractVectorField}
    metric::Metric
    d::Int  # dimensionality
    periodic_dimensions::Vector # label for which dimensions are periodic
end

Ring() = Ring(
    "ring",
    [0],
    [2π],
    [ConstantVectorField(1, 1)],
    PeriodicEuclidean([2π]), 1, [1]
)

apply_boundary_conditions!(x, ::Ring) = mod.(x, 2π), 1


# ---------------------------------------------------------------------------- #
#                                     PLANE                                    #
# ---------------------------------------------------------------------------- #
@with_repr struct Manifoldℝ² <: AbstractManifold
    name::String
    xmin::Vector
    xmax::Vector
    ψs::Vector{AbstractVectorField}
    metric::Metric
    d::Int  # dimensionality
    periodic_dimensions::Vector # label for which dimensions are periodic
end

Manifoldℝ²(m) = Manifoldℝ²(
    "ℝ²",
    [-m, -m],
    [m, m],
    [ConstantVectorField(2, 1), ConstantVectorField(2, 2)],
    Euclidean(),
    2, [0, 0]
)
ℝ² = Manifoldℝ²(100)

function apply_boundary_conditions!(x::Vector, m::Manifoldℝ²)
    vel_correction_factors = [1, 1]

    # non periodic dimension
    δ = 0.25  # padding around boundary to account for bump size
    if x[1] <= m.xmin[1] + δ
        x[1] = m.xmin[1] + δ
        vel_correction_factors[1] = 0
    elseif x[1] >= m.xmax[1] - δ
        x[1] = m.xmax[1] - δ
        vel_correction_factors[1] = 0
    end


    if x[2] <= m.xmin[2] + δ
        x[2] = m.xmin[2] + δ
        vel_correction_factors[2] = 0
    elseif x[2] >= m.xmax[2] - δ
        x[2] = m.xmax[2] - δ
        vel_correction_factors[2] = 0
    end

    return x, vel_correction_factors
end

function Base.rand(m::Manifoldℝ²; δ=0.25)
    # generate a random point on the unit sphere
    d = length(m.xmin)
    x = zeros(d)
    for i = 1:d
        x[i] = rand(m.xmin[i]+δ:0.001:m.xmax[i]-δ)  # padding because of boundary effect on neural activity
    end
    x
end



function upsample(x, n_steps)
    t = 1:length(x)
    itp = LinearInterpolation(t, x)
    new_t = range(1, length(x), length=n_steps)
    return itp.(new_t)
end



function make_space_filling_trajectory(M::Manifoldℝ², n_pts::Int, n_steps::Int)
    # sample a grid of points on the torus
    δ = 1
    x_min, x_max = M.xmin[1] + δ, M.xmax[1] - δ
    y_min, y_max = M.xmin[2] + δ, M.xmax[2] - δ
    x = range(x_min, x_max, length=n_pts)
    y = range(y_min, y_max, length=n_pts)

    # Create grid of coordinates
    points = Vector{Vector{Float64}}()
    
    # Create a snake-like pattern to ensure continuity
    for (i, xi) in enumerate(x)
        ys = i % 2 == 1 ? y : reverse(y)
        for yi in ys
            push!(points, [xi, yi])
        end
    end

    # Convert to matrix format (T×d)
    X = reduce(hcat, points)' |> Matrix
    println("X before upsampling: ", size(X))

    # Convert to matrix format (T×d)
    X = reduce(hcat, points)' |> Matrix
    X_up = zeros(n_steps, 2)
    X_up[:, 1] = upsample(X[:, 1], n_steps)
    X_up[:, 2] = upsample(X[:, 2], n_steps)
    println("X after upsampling: ", size(X_up), " with n_steps: ", n_steps)
    

    # Compute velocities between consecutive points
    V = similar(X_up)
    for i in 1:size(X_up,1)-1
        V[i,:] = (X_up[i+1,:] - X_up[i,:])/0.5
    end
    V[end,:] .= 0  # No velocity for last point
    
    return X_up, V
end


# ---------------------------------------------------------------------------- #
#                                   CYLINDER                                   #
# ---------------------------------------------------------------------------- #
@with_repr struct Cylinder <: AbstractManifold
    name::String
    xmin::Vector
    xmax::Vector
    ψs::Vector{AbstractVectorField}
    metric::Metric
    d::Int  # dimensionality
    periodic_dimensions::Vector # label for which dimensions are periodic
end

Cylinder(extent=5) = Cylinder(
    "Cylinder",
    [0, -extent],
    [2π, extent],
    [ConstantVectorField(2, 1), ConstantVectorField(2, 2)],
    PeriodicEuclidean([2π]),
    2, [1, 0]
)

C = Cylinder()

function apply_boundary_conditions!(x::Vector, m::Cylinder)
    
    vel_correction_facors = [1, 1]
    # non periodic dimension
    δ = 0.5  # padding around boundary to account for bump size
    if x[2] <= m.xmin[2] + δ
        x[2] = m.xmin[2] + δ
        vel_correction_facors[2] = 0
    elseif x[2] >= m.xmax[2] - δ
        x[2] = m.xmax[2] - δ
        vel_correction_facors[2] = 0
    end

    # periodic dimension
    x[1] = mod(x[1], 2π)
    return x, vel_correction_facors

end

function Base.rand(m::Cylinder; δ=0.25)
    # generate a random point on the unit sphere
    d = length(m.xmin)
    x = zeros(d)
    for i = 1:d
        _δ = d == 1 ? δ : 0.0
        x[i] = rand(m.xmin[i]+_δ:0.001:m.xmax[i]-_δ)  # padding because of boundary effect on neural activity
    end
    x
end

# ---------------------------------------------------------------------------- #
#                                     TORUS                                    #
# ---------------------------------------------------------------------------- #
@with_repr struct Torus <: AbstractManifold
    name::String
    xmin::Vector
    xmax::Vector
    ψs::Vector{AbstractVectorField}
    metric::Metric
    d::Int  # dimensionality
    periodic_dimensions::Vector # label for which dimensions are periodic
end

Torus() = Torus(
    "Torus",
    [0, 0],
    [2π, 2π],
    [ConstantVectorField(2, 1), ConstantVectorField(2, 2)],
    PeriodicEuclidean([2π, 2π]), 2, [1, 1]
)
T = Torus()

function apply_boundary_conditions!(x::Vector, ::Torus)
    return mod.(x, 2π), ones(length(x))
end


# ---------------------------------------------------------------------------- #
#                                    SPHERE                                    #
# ---------------------------------------------------------------------------- #
@with_repr struct Sphere <: AbstractManifold
    name::String
    xmin::Vector
    xmax::Vector
    ψs::Vector{AbstractVectorField}
    metric::Metric
    d::Int  # dimensionality
    periodic_dimensions::Vector # label for which dimensions are periodic
end
S² = Sphere(
    "S²",
    [-1, -1, -1],
    [1, 1, 1],
    [
        VectorField(sphere_ψx), 
        VectorField(sphere_ψy), 
        VectorField(sphere_ψz)
    ],
    SphericalDistance(),
    3, [0, 0, 0]
)

function apply_boundary_conditions!(x::Vector, ::Sphere)
    return x ./ norm(x), ones(length(x))
end

function Base.rand(::Sphere; δ=0)
    # generate a random point on the unit sphere
    x = rand(3) .- 0.5
    return x ./ norm(x)
end


# ---------------------------------------------------------------------------- #
#                                    MOBIUS                                    #
# ---------------------------------------------------------------------------- #
@with_repr struct Mobius <: AbstractManifold
    name::String
    xmin::Vector
    xmax::Vector
    ψs::Vector{AbstractVectorField}
    metric::Metric
    d::Int  # dimensionality
    periodic_dimensions::Vector # label for which dimensions are periodic
end

Mobius(extent=2) = Mobius(
    "Mobius",
    [-extent, 0],
    [extent, 2π],
    [VectorField(MB_ψ1), VectorField(MB_ψ2)],
    MobiusEuclidean(),
    2, [0, 1]
)


""" 
Correct the position vector x to ensure
that it's on the MB. If it's too much to the side
on the non periodic dimension, put it at the manifold's boundary,
if it's along the periodic dimension gets its position module 2π.
"""
function apply_boundary_conditions!(x::Vector, m::Mobius)
    vel_correction_factors = [1, 1]
    # non periodic dimension
    δ = 0.5  # padding around boundary to account for bump size
    if x[1] <= m.xmin[1] + δ
        x[1] = m.xmin[1] + δ
        vel_correction_factors[1] = 0
    elseif x[1] >= m.xmax[1] - δ
        x[1] = m.xmax[1] - δ
        vel_correction_factors[1] = 0
    end

    # periodic dimension
    if x[2] >= 2π
        x[2] = x[2] - 2π
        x[1] = -x[1]
        vel_correction_factors = [-1, 0]
    elseif x[2] <= 0
        x[2] = 2π + x[2]
        x[1] = -x[1]
        vel_correction_factors = [-1, 0]
    else
        vel_correction_factors = [1, 1] 
    end
    return x, vel_correction_factors
end

function Base.rand(m::Mobius; δ=0.25)
    # generate a random point on the unit sphere
    d = length(m.xmin)
    x = zeros(d)
    for i = 1:d
        _δ = d == 1 ? δ : 0.0
        x[i] = rand(m.xmin[i]+_δ:0.001:m.xmax[i]-_δ)  # padding because of boundary effect on neural activity
    end
    x
end


# ---------------------------------------------------------------------------- #
#                                 Klein Botttle                                #
# ---------------------------------------------------------------------------- #


@with_repr struct KleinBottle <: AbstractManifold
    name::String
    xmin::Vector
    xmax::Vector
    ψs::Vector{AbstractVectorField}
    metric::Metric
    d::Int  # dimensionality
    periodic_dimensions::Vector # label for which dimensions are periodic
end

KleinBottle() = KleinBottle(
    "Klein Bottle",
    [0, 0],
    [2π, 2π],
    [ConstantVectorField(2, 1), ConstantVectorField(2, 2)],
    KleinBottleEuclidean(),
    2, 
    [1, 1]
)


""" 
Correct the position vector x to ensure
that it's on the Klein Bottle. If it's too much to the side

"""
function apply_boundary_conditions!(x::Vector, m::KleinBottle)
    vel_correction_factors = [1, 1]
    # non periodic dimension
    # δ = 0.5  # padding around boundary to account for bump size
    # if x[1] <= m.xmin[1] + δ
    #     x[1] = m.xmin[1] + δ
    #     vel_correction_factors[1] = 0
    # elseif x[1] >= m.xmax[1] - δ
    #     x[1] = m.xmax[1] - δ
    #     vel_correction_factors[1] = 0
    # end

    # periodic dimension
    if x[2] >= 2π
        x[2] = x[2] - 2π
        x[1] = -x[1]
    elseif x[2] <= 0
        x[2] = 2π + x[2]
        x[1] = -x[1]
    end
    return x, [1, 1] 
end

function Base.rand(m::KleinBottle; δ=0)
    rand(0:0.001:2π, 2)
end

# ---------------------------------------------------------------------------- #
#                                 COVER SPACES                                 #
# ---------------------------------------------------------------------------- #

@with_repr struct CoverSpace
    M::AbstractManifold  # variable manifold (cover space)
    N::AbstractManifold  # neural manifold (covered)
    ρ::Function          # cover map | it's ρ (\rho) not 'p'
    ρⁱ::Function         # inverse of the cover map
    λs::Vector{Function}  # functions with scaling from points in M to points in N
end


CoverSpace(M, N, ρ, ρⁱ) = CoverSpace(M, N, ρ, ρⁱ, repeat([identity], length(M.xmin)))
CoverSpace(M, N) = CoverSpace(M, N, identity, identity, repeat([identity], length(M.xmin)))
CoverSpace(M) = CoverSpace(M, M)

identity(x) = x


end
