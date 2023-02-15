function line_maker(
    cantype;
    n::Int = 256,
    offset_size = 0.15,
    σ = :softrelu,
    α = 25,
    k = LocalGlobalKernel(α = 1, σ = 1)
)
    # neurons position and distance function
    n = (n,)  # number of neurons in the ring

    # neurons coordinates and metric
    M = Line()
    ξ(i::Int)::Vector = [lerp(i, n[1], M.xmin[1], M.xmax[1])]  # neurons coordinates function
    d = Euclidean()  # distance function

    # cover map
    cover = CoverSpace(M)

    # make network
    return if cantype == :single
        SingleCAN(
            "line",
            cover,
            n,
            ξ,
            d,
            k;
            σ = σ,
        )
    else
        CAN("line", 
        cover, n, ξ, d, k; 
        offset_size = offset_size, 
        σ = σ, 
        α = α
        ) 

    end
end
