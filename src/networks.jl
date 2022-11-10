# @info "Creating standard attractor networks"

# # ---------------------------------------------------------------------------- #
# #                                 RING ATTRACTOR                               #
# # ---------------------------------------------------------------------------- #
# @info "creating ring attractor"
# # neurons position and distance function
# n = (256,)  # number of neurons in the ring
# ξ_r(i::Int)::Vector = [lerp(i, n[1], 0.0, 2π)]  # neurons coordinates function
# d_r = PeriodicEuclidean([2π])  # distance function

# # construct network
# ring_attractor = CAN("ring", n, ξ_r, d_r, MexicanHatKernel(); offset_size = 0.25, α = 0.1);

# --------------------------------------------------------------75-------------- #
#                                TORUS ATTRACTOR                               #
# # ---------------------------------------------------------------------------- #
# @info "creating torus attractor"
# # neurons position and distance function
# n = (64, 64)
# function ξ_t(i::Int, j::Int)::Vector  # neurons coordinates function
#     n̂_i, n̂_j = Int(n[1] / 2), Int(n[2] / 2)
#     [lerp(i, n[1], -n̂_i, n̂_i), lerp(j, n[2], -n̂_j, n̂_j)]   # ∈ [-n/2, n/2] × [-n/2, n/2]
# end
# d_t = PeriodicEuclidean([n...])  # distance function over a torus manifold
# # connectivity kernel

# k_t = DiffOfExpKernel(; λ = 13.0)  # ? used for hexagonal grid pattern
# # k_t = LocalGlobalKernel(α = 0.04, σ = 20.0, β = 0.025)  # ? used for single bump pattern

# # construct network
# cover = CoverSpace(ℝ², T, (x, y) -> [mod(x - 32, 64), mod(y - 32, 64)])
# torus_attractor = CAN("torus", cover, n, ξ_t, d_t, k_t;)  # ? if using DiffOfExpKernel, change to α=0.10315

# # ---------------------------------------------------------------------------- #
# #                               MOBIUS ATTRACTOR                               #
# # ---------------------------------------------------------------------------- #
# @info "creating mobius attractor"
# n = (64, 64)
# ξ_m(i::Int, j::Int)::Vector = [lerp(i, n[1], 0.0, 2π), lerp(j, n[2], 0.0, 2π)]  # ∈ [0, 2π] × [0, 2π]

# d_m = MobiusEuclidean(2π)

# # connectivity kernel
# k_m = DiffOfExpKernel(; λ = 1.5)

# mobius_attractor = CAN("mobius", n, ξ_m, d_m, k_m; offset_size = 0.2, α = 0.1)


# ---------------------------------------------------------------------------- #
#                                    SPHERE                                    #
# ---------------------------------------------------------------------------- #
# import .Manifolds: sphere_embedding

# @info "Creating sphere attractor"
# n = (64, 64)
# function ξ_s(i::Int, j::Int)::Vector
#     [lerp(i, n[1], -π, π), lerp(j, n[2], -π / 2, π / 2)]
# end
# d_s = SphericalAngle()
# k_s = DiffOfExpKernel(; λ = 0.75)

# cover = CoverSpace(S², S², (x, y) -> [mod(x, 1), mod(y, 1)])
# sphere_attractor = CAN("sphere", cover, n, ξ_s, d_s, k_s; 
#     offset_size=[0.2, 0.2, 0.1, 0.1], 
#     φ=sphere_embedding)

