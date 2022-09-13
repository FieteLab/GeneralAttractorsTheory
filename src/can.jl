module Can
    import Base.Iterators: product as ×  # cartesian product
    import Distances: Metric, pairwise
    import StaticArrays: SVector, SA_F64, SMatrix
    using Term.Progress
    using Plots
    import LinearAlgebra: ⋅, I
    
    
    export AbstractCAN, CAN

    using ..Kernels: AbstractKernel


    # --------------------------- activation functions --------------------------- #
    relu(x) = max(0, x)

    activations::Dict{Symbol, Function} = Dict(
        :relu => relu,
        :tanh => tanh,
    )

    # --------------------------------- abstract --------------------------------- #
    abstract type AbstractCAN end


    # ---------------------------------------------------------------------------- #
    #                                      CAN                                     #
    # ---------------------------------------------------------------------------- #
    mutable struct CAN <: AbstractCAN
        name::String
        n::NTuple{N,Int} where N           # number of neurons in each dimension
        d::Int                             # number of dimensions
        I::Vector{Tuple}                   # index (i,j...) of each neuron in the lattice
        X::Matrix                          # N × n_neurons matrix with coordinates of each neuron in lattice
        Ws::Vector{Array}                  # connectivity matrices with lateral offsets | length N
        kernel::AbstractKernel             # connectivity kernel
        σ::Function                        # activation function
        A::Matrix{Float64}                 # K×D map projecting v∈Φ to manifold dimensions
    end

    Base.string(can::CAN) = "CAN (dim=$(length(can.n))) - n neurons: $(can.n)"
    Base.print(io::IO, can::CAN) = print(io, string(can))
    Base.show(io::IO, ::MIME"text/plain", can::CAN) = print(io, string(can))

    function CAN(
            name::String,
            n::NTuple{N,Int},
            ξ::Function,
            metric::Metric,
            kernel::AbstractKernel;
            σ::Union{Symbol, Function}=:relu,
            offset_size::Number=1.0,                        # magnitude of weights offset (distance)
            offsets::Union{Nothing, Matrix} = nothing,      # offset directions, rows Aᵢ of A
            α::Float64=1.0,                                   # scales A matrix
        ) where N

        d = length(n)

        # check that ξ has the right form
        nargs = first(methods(ξ)).nargs - 1
        @assert N == nargs  "ξ should accept $(N) arguments, accepts: $(nargs)"
    
        rtype = Base.return_types(ξ, NTuple{N, Int})[1]
        @assert rtype <: AbstractVector "ξ should return a Vector with neuron coordinates, not $rtype"
    
        # get the index of every neuron in the lattice | Array of size `n`
        lattice_idxs::AbstractArray{NTuple{N, Int}} = ×(map(_n -> 1:_n, n)...) |> collect |> vec
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
            D =  pairwise(metric, X .- offset_size .* θ, X)

            # get connectivity matrix with kernel
            push!(Ws, kernel.k.(D ))
        end

        # get connectivity function
        σ = σ isa Symbol ? activations[σ] : σ
        
        @debug "ready" n lattice_idxs eltype(lattice_idxs) X eltype(X) typeof(Ws) eltype(Ws)
        return CAN(name, n, d, lattice_idxs, X, Ws, kernel, σ, α .* Matrix(hcat(offsets...)'))
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



    function get_offsets(offsets::Union{Nothing, Matrix}, d::Int, n::Tuple)::Vector
        isnothing(offsets) && (offsets = Matrix(1.0I, d, d))
        D, K = size(offsets)
        @assert D == d "Offsets matrix `A` should have $(d) rows, not $D"

        v₀ = ones(K)
        

        basevec(i, x, d) = begin
            î = (Int ∘ ceil)(i / 2)
            bv = zeros(d) 
            bv[î] = x
            bv
        end
        
        offsets = map(
            i -> basevec(i, Aᵢ(offsets, i) ⋅ v₀, d), 1:2d
        ) |> collect
        @assert length(offsets)==2length(n) "Expected $(2length(n)) got $(length(offsets)) offsets"
        @assert length(offsets[1])==d

        @debug "Making CAN given offsets" size(offsets) typeof(offsets) v₀ Aᵢ d n offsets
        return offsets
    end

end