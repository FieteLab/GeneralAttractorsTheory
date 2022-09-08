module Can
    import Base.Iterators: product as ×  # cartesian product
    import Distances: Metric, pairwise
    import StaticArrays: SVector, SA_F64, SMatrix
    using Term.Progress
    using Plots
    
    

    export AbstractCAN, CAN, Kernel


    # ---------------------------------------------------------------------------- #
    #                                    KERNEL                                    #
    # ---------------------------------------------------------------------------- #

    """
        struct Kernel
            k::Function
        end

    A `Kernel` represents a function used to compute the connectivity strength
    between neurons x₁ -> x₂ based on their distance Δx. 
    It's entirely defined by a scalar function `k`.

    The struct is a convenience structure used to ensure that `k` has the appropriate
    signature and to help with multiple dispatch (e.g. for plotting).
    """
    struct Kernel
        k::Function
        function Kernel(k::Function)
            nargs = first(methods(k)).nargs - 1
            @assert nargs == 1 "Kernel functions can only have one argument, not $nargs"
    
            rtype = Base.return_types(k, (Number,))[1]
            @assert rtype isa Union{Any, Float64}  "Kernel function should return a scalar, not $rtype"
            
            new(k)
        end
    end
    
    # make `K` callable
    (K::Kernel)(x::Float64)::Float64 = K.k(x)
    (K::Kernel)(x::Vector)::Vector = K.k.(x)

    Base.string(K::Kernel) = "Kernel: $(K.k)"
    Base.print(io::IO, k::Kernel) = print(io, string(k))
    Base.show(io::IO, ::MIME"text/plain", k::Kernel) = print(io, string(k))



    # ---------------------------------------------------------------------------- #
    #                                      CAN                                     #
    # ---------------------------------------------------------------------------- #

    abstract type AbstractCAN end

    mutable struct CAN <: AbstractCAN
        n::NTuple{N,Int} where N           # number of neurons in each dimension
        I::Vector{Tuple}                   # index (i,j...) of each neuron in the lattice
        X::Matrix                          # N × n_neurons matrix with coordinates of each neuron in lattice
        Ws::Vector{Array}                  # connectivity matrices with lateral offsets | length N
        kernel::Kernel                     # connectivity kernel
    end

    Base.string(can::CAN) = "CAN (dim=$(length(can.n))) - n neurons: $(can.n)"
    Base.print(io::IO, can::CAN) = print(io, string(can))
    Base.show(io::IO, ::MIME"text/plain", can::CAN) = print(io, string(can))

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


    function CAN(
            n::NTuple{N,Int},
            ξ::Function,
            d::Metric,
            kernel::Kernel;
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
            # push!(Ws, kernel.k.(D ))
            push!(Ws, D)
        end
        
        @debug "ready" n lattice_idxs eltype(lattice_idxs) X eltype(X) typeof(Ws) eltype(Ws)
        return CAN(n, lattice_idxs, X, Ws, kernel)
    end


end