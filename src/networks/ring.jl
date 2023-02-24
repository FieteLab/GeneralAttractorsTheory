function ring_maker(
    cantype;
    n::Int = 256,
    offset_size = 0.15,
    σ = :softrelu,
    α = 25,
    k = LocalGlobalKernel(α = 1, σ = 1),
    cover_manifold = :default,
    kwargs...
)
    # neurons position and distance function
    n = (n,)  # number of neurons in the ring

    # neurons coordinates and metric
    ξ(i::Int)::Vector = [lerp(i, n[1], 0.0, 2π - 2π / n[1])]  # neurons coordinates function
    d = PeriodicEuclidean([2π])  # distance function

    # cover map
    if cover_manifold == :default
        cover = CoverSpace(Ring())
    else
        ρ(x::Number) = x
        ρ(x) = x[1]

        ρⁱ(x::Number) = mod(x, 2π)
        ρⁱ(x) = mod(x[1], 2π)

        cover = CoverSpace(
            Line(; d=Real(4π)),  # goes from -4π to 4π - 4x cover
            Ring(),
            ρ, ρⁱ,
        )
    end

    # make network
    return if cantype == :single
        SingleCAN(
            "ring",
            cover,
            n,
            ξ,
            d,
            k;
            σ = σ,
        )
    else
        CAN("ring", 
        cover, n, ξ, d, k; 
        offset_size = offset_size, 
        σ = σ, 
        α = α
        ) 

    end
end
