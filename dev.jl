using GeneralAttractors
using LinearAlgebra: I, ⋅  # import identity matrix and dot product
# using Tullio
using BenchmarkTools: @benchmark


τ = 10  # milliseconds
can = torus_attractor


function step!(S::Array{Float64}, Ṡ::AbstractVector{Float64},  v::Vector{Float64}, σ::Function)
    a = sum(S, dims=2)
    for i in 1:2can.d
        # get state and connection matrix.
        # Sᵢ::Vector{Float64} = view(S, :, i)
        Sᵢ::Vector{Float64} = S[:, i]  # use view to reduce allocations

        # use view to reduce allocations
        W::Matrix{Float64} = can.Ws[i]

        # get effect of all copies of neurons
        # for j in 1:2can.d
        #     Sⱼ = S[:, j]
        #     Ṡ .+= W * Sⱼ  # ! ALLOCATIONS HEAVY
        # end        
        Ṡ .+= W * a

        # get effect of bias and velocity input
        a_idx = (Int ∘ ceil)(i / 2)
        Aᵢ::Vector{Float64} = view(A, :, a_idx)
        Ṡ .+= b₀ + Aᵢ ⋅ v  # bias + velocity input

        # update S
        S[:, i] += σ.((Ṡ .- Sᵢ)./τ)

        # rest Ṡ
        Ṡ .*= 0.0
    end
    # return S
end

# ------------------------------------ run ----------------------------------- #
N = *(can.n...)
S = rand(N, 2can.d)

# allocate Ṡ to speed up
Ṡ = zeros(N)


# define Φ → can mapping matrix
K = 4  # dimensionality of variable space
A = Matrix(1.0I, K, can.d)

# get input velocity vector
v = rand(K)

# bias
b₀ = 0.0

step!(S, Ṡ, v, can.σ)
print("\n\n\n")
# @btime step!($S, $Ṡ, $v, $can.σ) 

@time for i in 1:100
    step!(S, Ṡ, v, can.σ)
end


# TODO improve performance
# TODO check validity
# TODO add visualization