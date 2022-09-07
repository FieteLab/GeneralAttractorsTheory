# ---------------------------------------------------------------------------- #
#                                 RING ATTRACTOR                               #
# ---------------------------------------------------------------------------- #

# neurons position and distance function
n = (256,)  # number of neurons in the ring
ξ_r(i::Int)::Vector = [(i-1)/(n[1]-1)*2π]  # neurons coordinates function
d_r = PeriodicEuclidean([2π])  # distance function

# connectivity kernel
α, σ = .33, .25
mexicanhat(t) = α * 2/(√(3σ)*π^(1/4))*(1 - (t/σ)^2) * exp(-t^2/(2σ))
k_r(x::Float64)::Float64 = mexicanhat(x)

# construct network
ring_attractor = CAN(n, ξ_r, d_r, Kernel(k_r); offset_strength=0.25);

# ---------------------------------------------------------------------------- #
#                                TORUS ATTRACTOR                               #
# ---------------------------------------------------------------------------- #
# neurons position and distance function
n = (48, 48)  
ξ_t(i::Int, j::Int)::Vector = begin # neurons coordinates function
    p_i, p_j = (i-1)/(n[1]-1), (j-1)/(n[2]-1) # ∈ [0, 1]
    [
        2π*(p_i), 
        2π*(p_j)
    ]  # ∈ [-π/2, π/2] × [-π/2, π/2]
end
d_t = PeriodicEuclidean([2π, 2π])  # distance function over a torus manifold

# connectivity kernel | params defined in fn to speedup
k_t(x::Float64)::Float64 = begin    
    a = 1.0
    λ = 13/2π
    β = 5/(λ^2)
    γ = 1.05 * β

    _x = abs(x)^2
    a * exp(-γ*_x) - exp(-β*_x)
end

# construct network
torus_attractor = CAN(n, ξ_t, d_t, Kernel(k_t); offset_strength=.5)