function ring_maker(
    cantype;
    n::Int = 256,
    offset_size = 0.15,
    σ = :softrelu,
    α = 25,
    k = LocalGlobalKernel(α = 1, σ = 1),
    kwargs...
)
    # neurons position and distance function
    n = (n,)  # number of neurons in the ring

    # neurons coordinates and metric
    ξ(i::Int)::Vector = [lerp(i, n[1], 0.0, 2π - 2π / n[1])]  # neurons coordinates function
    d = PeriodicEuclidean([2π])  # distance function

    # cover map
    cover = CoverSpace(Ring())

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
