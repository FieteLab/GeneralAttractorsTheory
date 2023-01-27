function plane_maker(
    cantype;
    n::Int = 48,
    k::AbstractKernel = LocalGlobalKernel(α = 2.5, σ = 5.0), 
    offset_size::Number = 0.2,
    α = 3.2,
    σ = :softrelu,
)

    # number of neurons
    n = (n, n) # number of neurons per dimension
    m_to_n(x) = x ./ 25
    n_to_m(x) = x .* 25

    ρ(x, y) = [x/25, y/25]
    ρ(v) = ρ(v...)
    ρⁱ(x, y) = [x*25, y*25]
    ρⁱ(w::Vector) = ρⁱ(w...)
    
    cover = CoverSpace(Manifoldℝ²(25), Manifoldℝ²(1), ρ, ρⁱ, [n_to_m, n_to_m])

    # define a function to get the coordinates of each neuron in the lattice
    function ξ(i::Int, j::Int)::Vector  # neurons coordinates function
        [lerp(i, n[1], -1, 1), lerp(j, n[2], -1, 1)]  
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


planecan = plane_maker(:can)
planecan_single = plane_maker(:single)