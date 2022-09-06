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
    # (K::Kernel)(x::Float64)::Float64 = K.k(x)
    (K::Kernel)(x::Vector)::Vector = K.k.(x)

    Base.string(k::Kernel) = "Kernel: $(K.k)"
    Base.print(io::IO, k::Kernel) = print(io, string(k))
    Base.show(io::IO, ::MIME"text/plain", k::Kernel) = print(io, string(k))



    # ---------------------------------------------------------------------------- #
    #                                      CAN                                     #
    # ---------------------------------------------------------------------------- #

    abstract type AbstractCAN end

    mutable struct CAN <: AbstractCAN
        n::NTuple{N,Int} where N           # number of neurons in each dimension
        I::Vector{Tuple{M, Int}} where M   # index (i,j...) of each neuron in the lattice
        X::Vector{Vector}                  # vector witht the position (vector) of each neuron in the neural lattice
        W::Array                           # connectivity matrix between all neurons in the network
    end

    function CAN(n::NTuple{N,Int}, ξ::Function, d::Metric, kernel::Kernel) where N
        # check that ξ has the right form
        nargs = first(methods(ξ)).nargs - 1
        @assert N == nargs  "ξ should accept $(N) arguments, accepts: $(nargs)"
    
        rtype = Base.return_types(ξ, NTuple{N, Int})[1]
        @assert rtype <: AbstractVector "ξ should return a Vector with neuron coordinates, not $rtype"
    
        # get the index of every neuron in the lattice | Array of size `n`
        lattice_idxs::AbstractArray{NTuple{N, Int}} = ×(map(_n -> 1:_n, n)...) |> collect |> vec
        @info "Got lattice" typeof(lattice_idxs) size(lattice_idxs)
    
        # get the coordinates of every neurons | Array of size `n`
        ξ̂(t::Tuple) = ξ(t...)
        X::Vector{Vector} = @time ξ̂.(lattice_idxs)
        @info "X" size(X) typeof(X) eltype(X) 
    
        # get connectivity matrix by applying the kernel function to the pairwise distance mtx
        D = @time  pairwise(d, hcat(X...))
        @info "got pairwise" size(D)
    
        W = @time kernel.k.(D)
        @info "Got connectivity" size(W)
    
        return CAN(n, lattice_idxs, X, W)
    end


end