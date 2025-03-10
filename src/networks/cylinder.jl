function cylinder_maker(
    cantype;
    n::Int = 48,
    k::AbstractKernel = LocalGlobalKernel(α = 2.5, σ = 25.0), 
    offset_size::Number = 0.2,
    α = 3.2,
    σ = :softrelu, 
    r2_extent = 25,
    cy_extent = 5,
    kwargs...
)

    # number of neurons
    n = (n, n) # number of neurons per dimension

    # ℝ² → Cy cover map.  
    ratio = r2_extent / cy_extent
    M = Manifoldℝ²(r2_extent)

    # define a function to map x ∈ -r2_extent, r2_extent to x ∈ -1,1
    r2_to_cy_scaling(x) = x / ratio
    # and the inverse
    cy_to_r2_scaling(x) = x * ratio

    """ ρ """
    ρ(x, y) = [mod(x, 2π), r2_to_cy_scaling(y)] # ρ: ℝ² → T
    ρ(v) = ρ(v...)


    """
        ρⁱ(x, y)

    Inverse of the cover map rho over the domain.
    Given a point (x,y) in N it gives a set of points (x̂, ŷ)
    in the cover space such that ρ(x̂, ŷ)=(x,y)
    """
    function ρⁱ(x, y; n = 40)
        ŷ = cy_to_r2_scaling(y)

        pts = zeros(2, n^2)
        for (c, i) in enumerate(-n/2:(n/2-1))
            x̂ = x + 2π * i
            pts[:, c*n] = [x̂, ŷ]
        end
        return pts
    end
    ρⁱ(w::Vector; n = 6) = ρⁱ(w...; n = n)

    cover = CoverSpace(M, Cylinder(cy_extent), ρ, ρⁱ, [identity, cy_to_r2_scaling])

    # define a function to get the coordinates of each neuron in the lattice
    function ξ(i::Int, j::Int)::Vector  # neurons coordinates function
        sep = 2π / n[1]
        [lerp(i, n[1], 0, 2π - sep), lerp(j, n[2], -cy_extent, cy_extent)]  
    end

    # define a distance metric
    d = PeriodicEuclidean([2π, Inf]) 

    # TODO offsets and stuff?


    return if cantype == :single 
        SingleCAN(
            "cylinder",
            cover,
            n,
            ξ,
            d,
            k;
            σ = σ,
        )
    else
        CAN(
        "cylinder",
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


# cylindercan = cylinder_maker(:can)
# cylindecan_single = cylinder_maker(:single)