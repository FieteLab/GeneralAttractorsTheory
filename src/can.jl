module Can
import Base.Iterators: product as ×  # cartesian product
import Distances: Metric, pairwise
import StaticArrays: SVector, SA_F64, SMatrix
using Term.Progress
using Plots
import LinearAlgebra: ⋅, I


export AbstractCAN, CAN

using ..Kernels: AbstractKernel
using..Manifolds: AbstractManifold

# ---------------------------------------------------------------------------- #
#                              DIFFERENTIAL FORMS                              #
# ---------------------------------------------------------------------------- #
"""
Differential forms are used to provide spatially varyi-ing inputs onto a network. 
"""

"""

Diferential one forms in d dimensions are defined by:
    ω = [f_1(x), f_2(x), ..., f_d(x)]
with x being the coordinates of a point on the manifold.
We want our one-forms to be aligned to basis one forms 
    dx_i = [0, ..., 1, ..., 0]
so we have 
    ω_i = [0, ..., f_i(x), ..., 0]

such that each one form is defined by `f` and `i`. 
"""
struct OneForm
    i::Int
    f::Function
end

"""
    (ω::OneForm)(x::Vector)

Evaluate a one form at a point `x`.
"""
function (ω::OneForm)(x::Vector)::Vector
    o = zeros(length(x))
    o[ω.i] = ω.f(x[ω.i])
    return o
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

activations::Dict{Symbol,Function} = Dict(:relu => relu, :tanh => tanh)


abstract type AbstractCAN end


"""
    mutable struct CAN <: AbstractCAN
        name::String
        M::AbstractManifold                # variable manifold
        n::NTuple{N,Int} where {N}         # number of neurons in each dimension
        d::Int                             # number of dimensions
        I::Vector{Tuple}                   # index (i,j...) of each neuron in the lattice
        X::Matrix                          # N × n_neurons matrix with coordinates of each neuron in lattice
        Ws::Vector{Array}                  # connectivity matrices with lateral offsets | length N
        kernel::AbstractKernel             # connectivity kernel
        σ::Function                        # activation function
        Ω::Vector{OneForm}                 # vector of `OneForm`s representing input measuring forms
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
    M::AbstractManifold                # variable manifold
    n::NTuple{N,Int} where {N}         # number of neurons in each dimension
    d::Int                             # number of dimensions
    I::Vector{Tuple}                   # index (i,j...) of each neuron in the lattice
    X::Matrix                          # N × n_neurons matrix with coordinates of each neuron in lattice
    Ws::Vector{Array}                  # connectivity matrices with lateral offsets | length N
    kernel::AbstractKernel             # connectivity kernel
    σ::Function                        # activation function
    Ω::Vector{OneForm}                 # vector of `OneForm`s representing input measuring forms
end

Base.string(can::CAN) = "CAN (dim=$(length(can.n))) - n neurons: $(can.n)"
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
    M::AbstractManifold,
    n::NTuple{N,Int},
    ξ::Function,
    metric::Metric,
    kernel::AbstractKernel;
    σ::Union{Symbol,Function} = :relu,
    offsets::Union{Nothing,Matrix} = nothing,         # offset directions, rows Aᵢ of A
) where {N}

    d = length(n)

    # check that ξ has the right form
    nargs = first(methods(ξ)).nargs - 1
    @assert N == nargs "ξ should accept $(N) arguments, accepts: $(nargs)"

    rtype = Base.return_types(ξ, NTuple{N,Int})[1]
    @assert rtype <: AbstractVector "ξ should return a Vector with neuron coordinates, not $rtype"

    # get the index of every neuron in the lattice | Array of size `n`
    lattice_idxs::AbstractArray{NTuple{N,Int}} = ×(map(_n -> 1:_n, n)...) |> collect |> vec
    @debug "Got lattice" typeof(lattice_idxs) size(lattice_idxs)

    # get the coordinates of every neurons | Array of size `n`
    ξ̂(t::Tuple) = ξ(t...)
    X::Matrix = hcat(ξ̂.(lattice_idxs)...)
    @debug "X" size(X) typeof(X) eltype(X)

    # get connectivity offset vectors
    offsets = get_offsets(offsets, d, n)

    # construct connectivity matrices
    Ws::Vector{Matrix} = []
    for θ in offsets
        # get pairwise offset connectivity
        D = pairwise(metric, X .- θ, X)

        # get connectivity matrix with kernel
        push!(Ws, kernel.k.(D))
    end

    # get connectivity function
    σ = σ isa Symbol ? activations[σ] : σ


    # construct one-forms
    Ω = OneForm[]

    # TODO construct Ω

    @debug "ready" n lattice_idxs eltype(lattice_idxs) X eltype(X) typeof(Ws) eltype(Ws)
    return CAN(name, M, n, d, lattice_idxs, X, Ws, kernel, σ, Ω)
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
function get_offsets(offsets::Union{Nothing,Matrix}, d::Int, n::Tuple)::Vector
    # initialize as an identiy matrix if note is provided
    isnothing(offsets) && (offsets = Matrix(1.0I, d, d))
    D, K = size(offsets)
    @assert D == d "Offsets matrix `A` should have $(d) rows, not $D"

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

    @debug "Making CAN given offsets" size(offsets) typeof(offsets) v₀ Aᵢ d n offsets
    return offsets
end

end
