using Plots
using GeneralAttractors
using SparseArrays
using Distances
using GeneralAttractors.Kernels
using GeneralAttractors.Manifolds
using GeneralAttractors: lerp

@info "creating torus attractor"
# neurons position and distance function
n = (64, 64)
function ξ_t(i::Int, j::Int)::Vector  # neurons coordinates function
    n̂_i, n̂_j = Int(n[1] / 2), Int(n[2] / 2)
    [lerp(i, n[1], -n̂_i, n̂_i), lerp(j, n[2], -n̂_j, n̂_j)]   # ∈ [-n/2, n/2] × [-n/2, n/2]
end
d_t = PeriodicEuclidean([n...])  # distance function over a torus manifold
# connectivity kernel

k_t = DiffOfExpKernel(; λ = 13.0)  # ? used for hexagonal grid pattern
# k_t = LocalGlobalKernel(α = 0.04, σ = 20.0, β = 0.025)  # ? used for single bump pattern

# construct network
cover = CoverSpace(ℝ², T, (x, y) -> [mod(x - 32, 64), mod(y - 32, 64)])
tor = CAN("torus", cover, n, ξ_t, d_t, k_t;)  # ? if using DiffOfExpKernel, change to α=0.10315



"""
Visualize a network's dynamics with each copy of the network beign independent
"""
Time = 200
can = tor
τ = 10

d = 2can.d
N = *(can.n...)
anim = Animation()

# prep matrices
S = [spzeros(Float64, N) for i in 1:d]
W = sparse.(map(x -> Float64.(x), can.Ws))
droptol!.(W, 0.001)

x̄ = range(0, maximum(can.X[1, :]), length = can.n[1])
ȳ = range(0, maximum(can.X[2, :]), length = can.n[2])

# simulate animate
for t in 1:Time
    t%10 == 0 && println("$t/$Time")

    # step
    for i in 1:d
        ṡ = W[i] * S[i] .+ 2
        S[i] += (can.σ.(ṡ)-S[i])/τ
    end

    t % 10 == 0 && begin
        plot(clims=(0, 1), levels=3)
        for i in 1:d
            offset = can.offsets[i] .* vec(maximum(can.X; dims = 2)) 
            x = offset[1] .+ x̄
            y = offset[2] .+ ȳ
            contourf!(x, y, reshape(S[i], can.n)', levels = 3)
        end
        frame(anim)
    end
    
end

gif(anim, "test.gif", fps = 10)