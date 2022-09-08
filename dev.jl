using GeneralAttractors


n = (48, 24)
function ξ_m(i::Int, j::Int)
    p_i, p_j = (i-1)/(n[1]-1), (j-1)/(n[2]-1) # ∈ [0, 1]
    [
        2π*(p_i), 
        p_j
    ]  # ∈ [0, 2π] × [0, 1]
end
d_m = MobiusEuclidean()

# connectivity kernel | params defined in fn to speedup
k_m(x::Float64)::Float64 = begin    
    a = 1.0
    λ = 5/2π
    β = 5/(λ^2)
    γ = 1.05 * β

    _x = abs(x)^2
    a * exp(-γ*_x) - exp(-β*_x)
end

mobius_attractor = @time CAN(n, ξ_m, d_m, Kernel(k_m); offset_strength=0.1)
show_connectivity(mobius_attractor)