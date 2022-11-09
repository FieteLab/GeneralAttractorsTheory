using Plots
using GeneralAttractors
using SparseArrays
using Distances
using GeneralAttractors.Kernels
using GeneralAttractors.Manifolds
using GeneralAttractors: lerp
using GeneralAttractors.Manifolds: sphere_embedding

# --------------------------------- make net --------------------------------- #
n = (64, 64)
function ξ_s(i::Int, j::Int)::Vector
    [lerp(i, n[1], -π, π), lerp(j, n[2], -π / 2, π / 2)]
end
d_s = SphericalAngle()
k_s = DiffOfExpKernel(; λ = 0.75)

cover = CoverSpace(S², S², (x, y) -> [x, y])
can = CAN("sphere", cover, n, ξ_s, d_s, k_s; 
        offset_size=[0.1, 0.1, 0.1, 0.1],
        φ=sphere_embedding
        )




# --------------------------------- simulate --------------------------------- #


"""
Visualize a network's dynamics with each copy of the network beign independent
"""
Time = 200
τ = 10
b = 1.0  # homogeneous input
η = 0.1  # noise std

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
        ṡ = W[i] * S[i] .+ b .+ (rand()-0.5) * η
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