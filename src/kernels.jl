module Kernels

export AbstractKernel, Kernel, MexicanHatKernel, DiffOfExpKernel, LocalGlobalKernel
export ConstantKernel
"""
    A `Kernel` represents a function used to compute the connectivity strength
        between neurons x₁ -> x₂ based on their distance Δx. 
        It's entirely defined by a scalar function `k`.
"""
abstract type AbstractKernel end

Base.string(K::AbstractKernel) = "Kernel: $(string(typeof(K)))"
Base.print(io::IO, k::AbstractKernel) = print(io, string(k))
Base.show(io::IO, ::MIME"text/plain", k::AbstractKernel) = print(io, string(k))


# make `K` callable
(K::AbstractKernel)(x::Float64)::Float64 = K.k(x)
(K::AbstractKernel)(x::Vector)::Vector = K.k.(x)


function check_kernel_function_signature(k::Function)
    nargs = first(methods(k)).nargs - 1
    @assert nargs == 1 "Kernel functions can only have one argument, not $nargs"

    rtype = Base.return_types(k, (Number,))[1]
    @assert rtype isa Union{Any,Float64} "Kernel function should return a scalar, not $rtype"
end


# ---------------------------------- Kernel ---------------------------------- #
"""
struct Kernel
    k::Function
end

Generic Kernel type, can be used with any function `k`.

A `Kernel` represents a function used to compute the connectivity strength
between neurons x₁ -> x₂ based on their distance Δx. 
It's entirely defined by a scalar function `k`.

The struct is a convenience structure used to ensure that `k` has the appropriate
signature and to help with multiple dispatch (e.g. for plotting).
"""
struct Kernel <: AbstractKernel
    k::Function
    function Kernel(k::Function)
        check_kernel_function_signature(k)
        new(k)
    end
end


# -------------------------------- mexican hat ------------------------------- #
"""
    MexicanHatKernel <: AbstractKernel
        α::Float64
        σ::Float64
        k::Function
    end

Kernel with a mexican hat function. α/σ define the shape.    
"""
struct MexicanHatKernel <: AbstractKernel
    α::Float64
    σ::Float64
    k::Function

    function MexicanHatKernel(; α = 0.33, σ = 0.25, β = 0)
        _k(t) = α * 2 / (√(3σ) * π^(1 / 4)) * (1 - (t / σ)^2) * exp(-t^2 / (2σ)) - β

        # get the peak location to ensure we can set it at 0
        peak = maximum(_k.(-4:.01:4))
        k(t) = _k(t) - peak
        new(α, σ, k)
    end
end


# -------------------------------- diff of exp ------------------------------- #

"""
    struct DiffOfExpKernel <: AbstractKernel
        a::Float64
        λ::Float64
        β::Float64
        γ::Float64
        k::Function

Kernel with a difference of exponential function.
Based on Burak & Fiete 2009.
"""
struct DiffOfExpKernel <: AbstractKernel
    a::Float64
    λ::Float64
    β::Float64
    γ::Float64
    k::Function

    function DiffOfExpKernel(;
        a::Float64 = 2.0,
        λ::Float64 = 5 / 2π,
        β::Float64 = 3 / (λ^2),
        γ::Float64 = 1.05 * β,
        δ::Float64 = 0.0,  # offset
    )

        _k(x) = begin
            _x = abs(x)^2
            a * exp(-γ * _x) - exp(-β * _x) + δ
        end

        peak = maximum(_k.(-4:.01:4))
        k(t) = _k(t) - peak
        new(a, λ, β, γ, k)
    end
end




# ------------------------------- local/global ------------------------------- #
"""
    struct LocalGlobalKernel <: AbstractKernel
        α::Float64   # strength of local excitation
        σ::Float64   # width of local excitation
        k:: Function

Kernel with local excitation and global inhibition. 
"""
struct LocalGlobalKernel <: AbstractKernel
    α::Float64   # strength of local excitation
    σ::Float64   # width of local excitation
    k::Function

    function LocalGlobalKernel(; α = 1.0, σ = 1.0)

        k(x) = α * exp(-(x^2 / 2σ * 2)) - α
        new(α, σ, k)
    end
end

"""
    ConstantKernel <: AbstractKernel
        σ::Float64  # half width of the kernel
        β_plus::Float64  # value of the kernel at x = 0
        β_minus::Float64  # value of the kernel at x >/< σ

Takes a value of β_plus at -σ < x < σ and β_minus otherwise.
"""
struct ConstantKernel <: AbstractKernel
    σ::Float64  # half width of the kernel
    β_plus::Float64  # value of the kernel at x = 0
    β_minus::Float64  # value of the kernel at x >/< σ
    k::Function

    function ConstantKernel(; σ = 1.0, β_plus = 0.0, β_minus = -1.0)
        k(x) = (x > σ || x < -σ) ? β_minus : β_plus
        new(σ, β_plus, β_minus, k)
    end
end

end # module Kernels



