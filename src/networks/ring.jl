function ring_maker(
    cantype;
    n::Int = 64,
    offset_size = 0.15,
    σ = :softrelu,
    α = 25,
    k = LocalGlobalKernel(α = 2.5, σ = 0.25)
)
    # neurons position and distance function
    n = (n,)  # number of neurons in the ring

    # neurons coordinates and metric
    ξ(i::Int)::Vector = [lerp(i, n[1], 0.0, 2π - 2π / n[1])]  # neurons coordinates function
    d = PeriodicEuclidean([2π])  # distance function

    # cover map
    cover = CoverSpace(Ring())

    # offsets and one forms
    # offsets = [
    #     p -> ring_ψ(p),
    #     p -> -ring_ψ(p)
    # ]

    # Ω = OneForm[OneForm(1, (x) -> ring_ψ(x)), OneForm(1, (x) -> -ring_ψ(x))]

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
        # offsets = offsets,
        # Ω = Ω,
        offset_size = offset_size, 
        σ = σ, 
        α = α
        ) 

    end
end

# ringcan = ring_maker(:can)
# ringcan_single = ring_maker(:single)