module Can
    import StaticArrays: SMatrix, SArray, SVector


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
    (K::Kernel)(x::Number) = K.k(x)

    Base.string(k::Kernel) = "Kernel: $(K.k)"
    Base.print(io::IO, k::Kernel) = print(io, string(k))
    Base.show(io::IO, ::MIME"text/plain", k::Kernel) = print(io, string(k))



    # ---------------------------------------------------------------------------- #
    #                                      CAN                                     #
    # ---------------------------------------------------------------------------- #

    abstract type AbstractCAN end

    mutable struct CAN <: AbstractCAN
        n::Tuple{Int}  # number of neurons in each dimension
        X::SVector     # vector witht the position (vector) of each neuron in the neural lattice
        W::SArray      # connectivity matrix between all neurons in the network
    end




end