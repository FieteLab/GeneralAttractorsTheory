function plane_maker(
    cantype;
    n::Int = 48,
    k::AbstractKernel = LocalGlobalKernel(α = 2.5, σ = 25.0), 
    offset_size::Number = 0.2,
    α = 3.2,
    σ = :softrelu, kwargs...
)

    # number of neurons
    n = (n, n) # number of neurons per dimension

    m_extent, n_extent = 25, 10 # size of the two manifolds
    ratio = m_extent / n_extent

    m_to_n(x) = x ./ ratio
    n_to_m(x) = x .* ratio

    ρ(x, y) = [x/ratio, y/ratio]
    ρ(v) = ρ(v...)
    ρⁱ(x, y) = [x*ratio, y*ratio]
    ρⁱ(w::Vector) = ρⁱ(w...)
    
    cover = CoverSpace(Manifoldℝ²(m_extent), Manifoldℝ²(n_extent), ρ, ρⁱ, [n_to_m, n_to_m])

    # define a function to get the coordinates of each neuron in the lattice
    function ξ(i::Int, j::Int)::Vector  # neurons coordinates function
        [lerp(i, n[1], -n_extent, n_extent), lerp(j, n[2], -n_extent, n_extent)]  
    end

    # define a distance metric
    d = Euclidean()
    
    return if cantype == :single 
        SingleCAN(
            "plane",
            cover,
            n,
            ξ,
            d,
            k;
            σ = σ,
        )
    else
        CAN(
        "plane",
        cover,
        n,
        ξ,
        d,
        k;
        offset_size = offset_size,
        σ = σ,
        α = α,
        # offsets = offsets,
        # Ω = Ω
    )
    end
end


# planecan = plane_maker(:can)
# planecan_single = plane_maker(:single)