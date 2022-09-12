# ---------------------------------------------------------------------------- #
#                                 RING ATTRACTOR                               #
# ---------------------------------------------------------------------------- #

# neurons position and distance function
n = (256,)  # number of neurons in the ring
ξ_r(i::Int)::Vector = [(i-1)/(n[1]-1)*2π]  # neurons coordinates function
d_r = PeriodicEuclidean([2π])  # distance function

# construct network
ring_attractor = CAN(n, ξ_r, d_r, MexicanHatKernel(); offset_strength=0.25);

# ---------------------------------------------------------------------------- #
#                                TORUS ATTRACTOR                               #
# ---------------------------------------------------------------------------- #
# neurons position and distance function
n = (64, 64)  
function ξ_t(i::Int, j::Int)::Vector  # neurons coordinates function
    p_i, p_j = (i-1)/(n[1]-1), (j-1)/(n[2]-1) # ∈ [0, 1]
    n̂_i, n̂_j = Int(n[1]/2), Int(n[2]/2)
    [
        -n̂_i*(1-p_i)+n̂_i*p_i, 
        -n̂_j*(1-p_j)+n̂_j*p_j, 
    ]  # ∈ [-n/2, n/2] × [-n/2, n/2]
end
d_t = PeriodicEuclidean([n...])  # distance function over a torus manifold
# connectivity kernel

k_t = DiffOfExpKernel(; λ = 13.0)

# construct network
torus_attractor = CAN(n, ξ_t, d_t, k_t; offset_strength=1)

# ---------------------------------------------------------------------------- #
#                               MOBIUS ATTRACTOR                               #
# ---------------------------------------------------------------------------- #

n = (48, 24)
function ξ_m(i::Int, j::Int)
    p_i, p_j = (i-1)/(n[1]-1), (j-1)/(n[2]-1) # ∈ [0, 1]
    [
        2π*(p_i), 
        p_j
    ]  # ∈ [0, 2π] × [0, 1]
end
d_m = MobiusEuclidean()

# connectivity kernel
k_m = DiffOfExpKernel(; λ = 5/2π)

mobius_attractor = CAN(n, ξ_m, d_m, k_m; offset_strength=0.1)
