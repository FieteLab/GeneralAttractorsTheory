module Networks
    import Distances: Metric, pairwise
   
    
    export AbstractNetwork, IntegratorNetwork, VelocityNetwork, CAN

    using ..Kernels: AbstractKernel


    # --------------------------- activation functions --------------------------- #
    relu(x) = max(0, x)

    activations::Dict{Symbol, Function} = Dict(
        :relu => relu,
        :tanh => tanh,
    )


    # ---------------------------------------------------------------------------- #
    #                            INTEGRATOR NETWORK (G)                            #
    # ---------------------------------------------------------------------------- #
    abstract type AbstractNetwork end

    Base.print(io::IO, net::AbstractNetwork) = print(io, string(net))
    Base.show(io::IO, ::MIME"text/plain", net::AbstractNetwork) = print(io, string(net))


    struct IntegratorNetwork <: AbstractNetwork
        n::NTuple{N,Int} where N           # number of neurons in each dimension
        d::Int                             # number of dimensions
        I::Vector{Tuple}                   # index (i,j...) of each neuron in the lattice
        X::Matrix                          # N × n_neurons matrix with coordinates of each neuron in lattice
        W::Vector{Array}                   # connectivity matrix (n × n)
        kernel::AbstractKernel             # connectivity kernel
        σ::Function                        # activation function
        τ::Float64                         # network time constant
    end

    Base.string(gnet::IntegratorNetwork) = "IntegratorNetwork (dim=$(gnet.d)) - n neurons: $(gnet.n)"



    """
        IntegratorNetwork(
            n::NTuple{N,Int},
            ξ::Function,
            metric::Metric,
            kernel::AbstractKernel;
            σ::Union{Symbol, Function}=:relu,
        ) where N

    Integrator network constructor. 
    To construct the network, get a matrix of pairwise distances `D` based on each 
    neuron's location in the lattice (given by `ξ` and `metric`). Then use a `kernel`
    to compute the connectivity strength based on distance. The `metric` function imposes
    a topology on the network connectivity.

    ## Arguments
        - n:            number of neurons along each dimension of the neural lattice.
        - ξ:            coordinates function. Gets the coordinates of each neuron in the lattice based on its index.
        - metric:       distance function (metric). Computes distance between neurons based on lattice topology.
        - kernel:       kernel type for weights matrix construction
        - σ:            non-linearity function
        - τ:            network time constant (milliseconds)
    """
    function IntegratorNetwork(
            n::NTuple{N,Int},
            ξ::Function,
            metric::Metric,
            kernel::AbstractKernel;
            σ::Union{Symbol, Function}=:relu,
            τ::Float64  = 10.0 
        ) where N

        # ---------------------------------- checks ---------------------------------- #
        d = length(n)

        # check that ξ has the right form
        nargs = first(methods(ξ)).nargs - 1
        @assert N == nargs  "ξ should accept $(N) arguments, accepts: $(nargs)"
    
        rtype = Base.return_types(ξ, NTuple{N, Int})[1]
        @assert rtype <: AbstractVector "ξ should return a Vector with neuron coordinates, not $rtype"
    
        # get the index of every neuron in the lattice | Array of size `n`
        lattice_idxs::AbstractArray{NTuple{N, Int}} = ×(map(_n -> 1:_n, n)...) |> collect |> vec
        @debug "Got lattice" typeof(lattice_idxs) size(lattice_idxs)
    
        # --------------------------------- construc --------------------------------- #
        # get the coordinates of every neurons | Array of size `n`
        ξ̂(t::Tuple) = ξ(t...)
        X::Matrix = hcat(ξ̂.(lattice_idxs)...)
        @debug "X" size(X) typeof(X) eltype(X) 
    
        
        # construct connectivity matrices (first get distance then pass through kernel)
        # get pairwise offset connectivity
        D = pairwise(metric, X, X)
        W = kernel.k.(D )
    
        # get connectivity function
        σ = σ isa Symbol ? activations[σ] : σ
        
        @debug "ready" n lattice_idxs eltype(lattice_idxs) X eltype(X) typeof(Ws) eltype(Ws)
        return IntegratorNetwork(n, d, lattice_idxs, X, W, kernel, σ, τ)
    end

    # ---------------------------------------------------------------------------- #
    #                             VELOCITY NETWORK (H)                             #
    # ---------------------------------------------------------------------------- #
    """
        A velocity input network of a CAN.

        An H net is composed of nᵢ neurons (neurons in the i-th dimension of a G network)
        indexed by θᵢ coordinates. It provides an input of the v⃗ᵢ velocity component to a
        G network. It receives inputs from G through a B matrix and provides input to it
        through an A matrix. 
        It doesn't have recurrent connectivity. 
        It has a function ϕ that convert v⃗ᵢ into a vector added to the network dynamics.

        H networks work by taking the G network state on the i-th dimension, shifting it
        and sending it back to G to cause a shift in the dynamics (with scaling due to speed).
        For each i-th dimension there's a positive and a negative shift network. 
    """
    struct VelocityNetwork <: AbstractNetwork
        m::Int          # number of neurons
        i::Int          # index of dimension (for multi-dimensional CAN's)
        sign::Int       # 1 if it's a positive shift network or 0 otherwise.

        A::Matrix       # n × m (m = # neurons in H, m = tot # neurons in G). H → G projections
        P::Matrix       # n × m. Projection matrix. P = Bᵀ
        S::Matrix       # m × m. Shift operator matrix on the state of H.
        B::Matrix       # m × n. G → H projection

        τ::Float64      # network time constant
    end

    """
        VelocityNetwork(G::IntegratorNetwork, i::Int, sign::Int)

    Network constructor. 
    Construct an H^±ᵢ network from a G net.    
    """
    function VelocityNetwork(G::IntegratorNetwork, i::Int, sign::Int, τ::Float64)
        m = G.n[i]                # neurons in H net
        n = reduce(*, G.n)        # tot neurons in G net

        # construct B matrix
        B =

        # construct A matrix
        P = B' |> collect
        S = 
        A = P*S

        return VelocityNetwork(m, i, sign, A, P, S, B, τ)
    end
    
    


    # ---------------------------------------------------------------------------- #
    #                                      CAN                                     #
    # ---------------------------------------------------------------------------- #
    
    abstract type AbstractCAN <: AbstractNetwork end
    
    """
        A Continous Attractor Network is made up of one integrator network G and 
        multiple velocity input networks Hs.
    """
    mutable struct CAN <: AbstractCAN
        name::String
        G::IntegratorNetwork
        Hs::Vector{VelocityNetwork}
    end

    Base.string(can::CAN) = "CAN $(can.name) | n. neurons: $(can.G.n)"

    function CAN(
        name::String,
        n::NTuple{N,Int},
        ξ::Function,
        metric::Metric,
        kernel::AbstractKernel;
        σ::Union{Symbol, Function}=:relu,
        τ_G = 10.0, 
    )

        G = IntegratorNetwork(n, ξ, metric, kernel; σ=σ, τ=τ_G)
    end

end