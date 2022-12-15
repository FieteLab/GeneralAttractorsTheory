module Can
import Base.Iterators: product as ×  # cartesian product
import Distances: Metric, pairwise
import LinearAlgebra: ⋅, I, norm
import ForwardDiff: jacobian


export AbstractCAN, CAN, offset_for_visual

import ..GeneralAttractors: by_column, softrelu
using ..Kernels: AbstractKernel
using ..ManifoldUtils: CoverSpace, area_deformation

# ---------------------------------------------------------------------------- #
#                              DIFFERENTIAL FORMS                              #
# ---------------------------------------------------------------------------- #
"""
Differential forms are used to provide spatially varyi-ing inputs onto a network. 
"""

"""

Diferential one forms in d dimensions over the neuronal lattice.

One forms are defined by:
    ω = [f_1(x), f_2(x), ..., f_d(x)]
with x being the coordinates of a point on the manifold.

Constant  one-forms are assumed to be aligned to basis one forms 
    dx_i = [0, ..., 1, ..., 0]
so we have 
    ω_i = [0, ..., f_i(x), ..., 0]
such that each one form is defined by `f` and `i`. 

More generally, `f` is a function of points on the lattice giving a one form
in coordinates representation and `i` is ignored.
"""
struct OneForm
    i::Int
    f::Function
end

Base.string(ω::OneForm) = "OneForm($(ω.i), $(ω.f))"
Base.print(io::IO, ω::OneForm) = print(io, string(ω))
Base.show(io::IO, ::MIME"text/plain", ω::OneForm) = print(io, string(ω))


"""
    (ω::OneForm)(x::Vector)

Evaluate a one form at a point `x`.
"""
function (ω::OneForm)(x::Vector)::Vector
    nargs = first(methods(ω.f)).nargs - 1

    if nargs == 1 && length(x) > 1
        o = zeros(length(x))
        o[ω.i] = ω.f(x[ω.i])
        return o
    elseif length(x) == 1
        # one dimensional CAN
        return [ω.f(x...)]
    else
        return ω.f(x...)
    end
end

"""
    (ω::OneForm)(x::Vector, v::Vector)::number

Evaluate ωₚ(v).
"""
function (ω::OneForm)(x::Vector, v::Vector)::Number
    return v ⋅ ω(x)
end

# ---------------------------------------------------------------------------- #
#                                      CAN                                     #
# ---------------------------------------------------------------------------- #

# --------------------------- activation functions --------------------------- #
relu(x) = max(0, x)

function sigmoid(x, scale = 1, steepness = 1, switch_point = 0)
    return scale / (1 + exp(-steepness * (x - switch_point)))
end


activations::Dict{Symbol,Function} = Dict(
    :relu => relu,
    :tanh => x -> (tanh(x) + 1),
    :softrelu => x -> softrelu(x),
    :step => x -> x >= 0 ? 4 : 0,
    :sigmoid => x -> sigmoid(x, 8, 2, 0),
)


abstract type AbstractCAN end


"""
    mutable struct CAN <: AbstractCAN
        name::String
        C::CoverSpace
        n::NTuple{N,Int} where {N}         # number of neurons in each dimension
        d::Int                             # number of dimensions
        I::Vector{Tuple}                   # index (i,j...) of each neuron in the lattice
        X::Matrix                          # N × n_neurons matrix with coordinates of each neuron in lattice
        Ws::Vector{Array}                  # connectivity matrices with lateral offsets | length N
        kernel::AbstractKernel             # connectivity kernel
        σ::Function                        # activation function
        Ω::Vector{OneForm}                 # vector of `OneForm`s representing input measuring forms
        offsets::Vector                    # vector of offset directions
    end

A Continous Attractor Network. 

Each CAN is a d-dimensional integrator network with nᵢ neurons in each dimension. Neurons are indexed
and are assigned coordinates in the neural lattice (`X`) given a coordinate function `ξ` that acts on 
the neurons' indices. 
Each CAN is comprised of 2d "copies" of a neural population, each with its connectivity matrix Wᵢ. 
The `Ws` are computed based on offset pairwise neuron coordinates. For each copy the offeset is in a direction
specified by `offsets` vectors (generally ± each basis direction) and scaled by an `offset_size` magnitude. 
The velocity inputs onto each copy are scaled by ωᵢ(v) with `ωᵢ` being a `OneForm` aligned to the i-th
offset direction but (potentially) varying in magnitude over the variable manifold. 
"""
mutable struct CAN <: AbstractCAN
    name::String
    C::CoverSpace
    n::NTuple{N,Int} where {N}         # number of neurons in each dimension
    d::Int                             # number of dimensions
    I::Vector{Tuple}                   # index (i,j...) of each neuron in the lattice
    X::Matrix                          # N × n_neurons matrix with coordinates of each neuron in lattice
    Ws::Vector{Array}                  # connectivity matrices with lateral offsets | length N
    kernel::AbstractKernel             # connectivity kernel
    σ::Function                        # activation function
    Ω::Vector{OneForm}                 # vector of `OneForm`s representing input measuring forms
    offsets::Vector
    offset_size::Any
    metric::Metric
    α::Number                           # scaling factor to make bump speed the same as input speed
end

Base.string(can::CAN) = "CAN $(can.name) (dim=$(length(can.n))) - n neurons: $(can.n)"
Base.print(io::IO, can::CAN) = print(io, string(can))
Base.show(io::IO, ::MIME"text/plain", can::CAN) = print(io, string(can))

"""
    function CAN(
        name::String,
        n::NTuple{N,Int},
        ξ::Function,
        metric::Metric,
        kernel::AbstractKernel;
        σ::Union{Symbol,Function} = :relu,
        offsets::Union{Nothing,Matrix} = nothing,         # offset directions, rows Aᵢ of A
        α::Float64 = 1.0,                                 # scales A matrix
    ) where {N}

Construct a d-dimensional network. 
"""
function CAN(
    name::String,
    C::CoverSpace,
    n::NTuple{N,Int},
    ξ::Function,
    args...;
    kwargs...,
) where {N}
    # check that ξ has the right form
    nargs = first(methods(ξ)).nargs - 1
    @assert N == nargs "ξ should accept $(N) arguments, accepts: $(nargs)"

    rtype = Base.return_types(ξ, NTuple{N,Int})[1]
    @assert rtype <: AbstractVector "ξ should return a Vector with neuron coordinates, not $rtype"

    # get the index of every neuron in the lattice | vector of length N `n` with elements of length d
    lattice_idxs::AbstractArray{NTuple{N,Int}} = ×(map(_n -> 1:_n, n)...) |> collect |> vec
    @debug "Got lattice" typeof(lattice_idxs) size(lattice_idxs)

    # get the coordinates of every neurons
    ξ̂(t::Tuple) = ξ(t...)
    X::Matrix = hcat(ξ̂.(lattice_idxs)...)  # matrix, size d × N
    @debug "X" size(X) typeof(X) eltype(X)

    # finalize
    return CAN(name, C, n, lattice_idxs, X, args...; kwargs...)
end



function CAN(
    name::String,
    C::CoverSpace,
    n::NTuple{N,Int},
    lattice_idxs::Vector,
    X::Matrix,
    metric::Metric,
    kernel::AbstractKernel;
    σ::Union{Symbol,Function} = :relu,
    Ω::Union{Nothing,Vector{OneForm}} = nothing,      # one forms for input velocity
    offsets::Union{Nothing,Vector} = nothing,         # offset directions, rows Aᵢ of A
    offset_size::Number = 1.0,
    φ::Union{Function,Nothing} = nothing,          # an embedding function if the distance over M is computed on an embedding of M in ℝᵐ
    α::Number = 1,
) where {N}
    d = length(n)

    # get connectivity offset vectors
    offsets::Vector{AbstractWeightOffset} = get_offsets(offsets, d, n)
    @debug "Got offsets" offsets eltype(offsets) offset_size

    # construct connectivity matrices
    Ws::Vector{Matrix} = []
    for offset in offsets
        # get pairwise offset connectivity
        D::Matrix = get_pairwise_distance_with_offset(offset, X, metric, offset_size, φ)

        # get connectivity matrix with kernel
        W = kernel.k.(D)

        # store connectivity matrix
        push!(Ws, W)
    end

    # get connectivity function
    σ = σ isa Symbol ? activations[σ] : σ

    # construct one-forms
    Ω = get_one_forms(Ω, offsets)

    @debug "ready" n lattice_idxs eltype(lattice_idxs) X eltype(X) typeof(Ws) eltype(Ws)
    return CAN(
        name,
        C,
        n,
        d,
        lattice_idxs,
        X,
        Ws,
        kernel,
        σ,
        Ω,
        offsets,
        offset_size,
        metric,
        α,
    )
end


# ---------------------------------------------------------------------------- #
#                                   ONE FORMS                                  #
# ---------------------------------------------------------------------------- #

"""
    get_one_forms(Nothing, offsets::Vector)::Vector{OneForm}

Construct default One Forms (dxᵢ).
"""
function get_one_forms(::Nothing, offsets::Vector)::Vector{OneForm}
    Ω = OneForm[]
    for (i, v) in enumerate(offsets)
        î = (Int ∘ ceil)(i / 2)
        ω = OneForm(î, x -> v.θ[î])
        push!(Ω, ω)
    end
    Ω
end

"""
validate one forms passed at CAN creation. 
"""
function get_one_forms(Ω::Vector{OneForm}, offsets)::Vector{OneForm}
    @assert length(Ω) == length(offsets)
    return Ω
end


# ---------------------------------------------------------------------------- #
#                             Offsets &  Distances                             #
# ---------------------------------------------------------------------------- #

abstract type AbstractWeightOffset end

"""
Compute pairwise distance between neurons after applying weights offsets.
"""
function get_pairwise_distance_with_offset end


"""
Get a vector shift representation for visualization (e.g. to offset connectivity matrices)
"""
function offset_for_visual end

# ----------------------------- constant offsets ----------------------------- #

"""
Offsets defined on the neuronal lattice, causing a constant shift of the coordinates
on the lattice.
"""
struct ConstantOffset <: AbstractWeightOffset
    θ::Vector
end

(o::ConstantOffset)(args...) = o.θ

""" normalize shift for visualizations """
function offset_for_visual(off::ConstantOffset)
    o = off.θ
    return o ./ (o .+ 0.01) .* sign.(o)
end


"""
    Aᵢ(A::Matrix, i::Int)

Get the correct row of `A` for copy `i`
of the data. This involves dealing with the 
different sizes of A and # of copies and the need
to invert the sign.
"""
function Aᵢ(A::Matrix, i::Int)::Vector
    î = (Int ∘ ceil)(i / 2)
    return if i % 2 == 0
        -A[:, î]
    else
        A[:, î]
    end
end


"""
    get_offsets(offsets::Union{Nothing,Matrix}, d::Int, n::Tuple)::vector

Construct a matrix representing the offsets in the connectivity matrix. 
The same offsets are used to construct the matrix `A` scaling velocity
inputs onto each population copu. 
"""
function get_offsets(::Nothing, d::Int, n::Tuple)::Vector{ConstantOffset}
    # initialize as an identiy matrix if note is provided
    offsets = Matrix(1.0I, d, d)
    _, K = size(offsets)
    v₀ = ones(K)

    """ basis vector for the i-th copy """
    basevec(i, x, d) = begin
        î = (Int ∘ ceil)(i / 2)
        bv = zeros(d)
        bv[î] = x
        bv
    end

    offsets = map(i -> basevec(i, Aᵢ(offsets, i) ⋅ v₀, d), 1:2d) |> collect
    @assert length(offsets) == 2length(n) "Expected $(2length(n)) got $(length(offsets)) offsets"
    @assert length(offsets[1]) == d
    return ConstantOffset.(offsets)
end

function get_offsets(offsets::Vector{Number}, d::Int, ::Tuple)::Vector{ConstantOffset}
    @assert all(length.(offsets) .== d)
    return ConstantOffset.(offsets)
end




function get_pairwise_distance_with_offset(
    offset::ConstantOffset,
    X::Matrix,
    metric::Metric,
    offset_size::Number,
    args...,
)::Matrix
    return pairwise(metric, X .- (offset_size .* offset.θ), X)
end


# ---------------------------- field based offsets --------------------------- #

struct FieldOffset <: AbstractWeightOffset
    displacement::Any  # used in-place of the real offsets e.g. to display connectivity mtxs during visualization
    ψ::Function
end
(o::FieldOffset)(p) = o.ψ(p)
(o::FieldOffset)(args...) = o(vec(args))


""" normalize shift for visualizations """
offset_for_visual(off::FieldOffset) = off.displacement


"""
Construct offsets given a list of vector field functions
"""
function get_offsets(offsets::Vector{Function}, ::Int, n::Tuple)::Vector{FieldOffset}
    @assert length(offsets) >= 2length(n) "Got $(length(offsets)) offset fields, expected at least: $(2length(n))"

    # get "base" offsets and use them as displacements FOR PLOTTING ONLY
    k = length(offsets)
    displacements = map(i -> 1.5 .* [cos(i * 2π / k), sin(i * 2π / k)], 1:k)

    return FieldOffset.(displacements, offsets)
end


"""
Get pairwise neurons distance given vector field offsets applied directly to the neurons lattice
"""
function get_pairwise_distance_with_offset(
    offset::FieldOffset,
    X::Matrix,
    metric::Metric,
    offset_size::Number,
    ::Nothing,
)::Matrix
    # δx = by_column(offset.ψ, X)

    # n = norm(δx)

    δx = hcat(
        map(
            x -> offset.ψ(x),
            eachcol(X)
        )...
    )

    pairwise(metric, X .- (offset_size .* δx), X)
    # X̂ = by_column(x -> x .- offset.ψ(x), X)
    # pairwise(metric, X̂, X)
end


"""
---
    get_pairwise_distance_with_offset(
        offset::FieldOffset,
        X::Matrix,
        metric::Metric,
        offset_size::Number,
        φ::Function,
    )::Matrix

Offset vectors in the pairwise distance are computed on an embedding of
the base manifold. Prior to computing the distance, neurons coordinates are 
offset according to a vector given by a vector field defined over the embedded manifold. 
"""
function get_pairwise_distance_with_offset(
    offset::FieldOffset,
    X::Matrix,
    metric::Metric,
    offset_size::Number,
    φ::Function,
)::Matrix
    # get coordinates in embedding space
    Y = by_column(φ, X)  # matrix q × n | q: embedding dimension, n: number of points

    # get vectors in embedding space
    V = by_column(offset.ψ, Y)

    # get vectors in manifold domain

    """ 
    Get vector in the manifold domain. 
    x∈X, v∈V. 

    w = Jᵀv ∈ ℝᵈ is a vector on the lattice of neurons
    """
    φᵀ(x::Vector, v::Vector) = jacobian(φ, x)' * v

    N = size(X, 2)
    W = hcat(map(i -> φᵀ(X[:, i], V[:, i]), 1:N)...)

    # apply shifts to neurons coordinates in the lattice & get distance
    pairwise(metric, X .- (offset_size .* W), X)
end




end
