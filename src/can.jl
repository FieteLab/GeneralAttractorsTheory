module Can
    import Base.Iterators: product as ×  # cartesian product
    import Distances: Metric, pairwise
    import StaticArrays: SVector, SA_F64, SMatrix
    using Term.Progress
    using Plots
    
    
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
        n::NTuple{N,Int} where N           # number of neurons in each dimension
        d::Int                             # number of dimensions
        I::Vector{Tuple}                   # index (i,j...) of each neuron in the lattice
        X::Matrix                          # N × n_neurons matrix with coordinates of each neuron in lattice
        Ws::Vector{Array}                  # connectivity matrices with lateral offsets | length N
        kernel::AbstractKernel             # connectivity kernel
        σ::Function                        # activation function
        offset_directions::Vector{Vector}  # offset direction for each copy of the neurons
    end

    Base.string(can::CAN) = "CAN (dim=$(length(can.n))) - n neurons: $(can.n)"
    Base.print(io::IO, can::CAN) = print(io, string(can))
    Base.show(io::IO, ::MIME"text/plain", can::CAN) = print(io, string(can))

    function CAN(
            n::NTuple{N,Int},
            ξ::Function,
            d::Metric,
            kernel::AbstractKernel;
            σ::Union{Symbol, Function}=:relu,
            offset_strength::Number=1.0
        ) where N

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
    
        # construct connectivity matrices with on vector offets
        orientations = make_orientations_table(n)
        Ws::Vector{Matrix} = []
        for θ in orientations
            # get pairwise offset connectivity
            D =  pairwise(d, X .+ offset_strength .* θ, X)

            # get connectivity matrix with kernel
            push!(Ws, kernel.k.(D ))
        end

        # get connectivity function
        σ = σ isa Symbol ? activations[σ] : σ
        
        @debug "ready" n lattice_idxs eltype(lattice_idxs) X eltype(X) typeof(Ws) eltype(Ws)
        return CAN(n, length(n), lattice_idxs, X, Ws, kernel, σ, orientations)
    end


    """
        make_orientations_table(::NTuple{N, Int})::Vector{Int} where N

    Given an N dimensional connectivity matrix, construct valid combinations
    of positive/negative basis vectors. 

    In ℝ, possible combinations are [-1], [1].
    In ℝ²: [-1, 0], [1, 0], [0, -1], [0, 1]
    In ℝ³: [1, 0, 0], [-1, 0, 0], [0, -1, 0], [0, 1, 0]...

    Avoid diagonal ([1, 1]) and 0 vectors. 

    These are used as offset vectors to create asymmetric connectivity matrix
    for attractors with drifting dynamics.
    """
    function make_orientations_table(::NTuple{N, Int})::Vector{Vector} where N
        v = [-1, 0, 1]  # possible values for each base vec of connection mtx
        θ = ×(repeat([v], N)...) |> collect  # all possible combinations
        θ = map(x -> [x...], vec(θ))  # turn into a vector of vectors

        # keep only elements of the form: [0, 1], [-1, 0]... ∈ ℝ² | [1, 1, 0], [-1, 0, 1] ∈ ℝ²....
        filter!(  
            x -> abs(sum(x)) == 1.0, 
            θ
        )
    end
end