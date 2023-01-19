
function torus_maker(cantype; 
        m::Int=64,
        k_t::AbstractKernel = LocalGlobalKernel(α = 2.5, σ = 5.0), 
        offset_size::Number = 0.2,
    )
    # number of neurons
    n = (m, m) # number of neurons per dimension

    # ℝ² → T cover map.
    """ ρ """
    ρ(x, y) = [mod(x, 2π), mod(y, 2π)] # ρ: ℝ² → T
    ρ(v) = ρ(v...)

    """
        ρⁱ(x, y)

    Inverse of the cover map rho over the domain.
    Given a point (x,y) in N it gives a set of points (x̂, ŷ)
    in the cover space such that ρ(x̂, ŷ)=(x,y)
    """
    function ρⁱ(x, y; n = 40)
        pts = zeros(2, n^2)
        for (c, i) in enumerate(-n/2:(n/2-1)), (k, j) in enumerate(-n/2:(n/2-1))
            x̂ = x + 2π * i
            ŷ = y + 2π * j
            pts[:, (c-1)*n+k] = [x̂, ŷ]
        end
        return pts
    end
    ρⁱ(w::Vector; n = 6) = ρⁱ(w...; n = n)

    cover = CoverSpace(Manifoldℝ²(100), Torus(), ρ, ρⁱ)

    # define a function to get the coordinates of each neuron in the lattice
    function ξ_t(i::Int, j::Int)::Vector  # neurons coordinates function
        sep = 2π / n[1]
        [lerp(i, n[1], 0, 2π - sep), lerp(j, n[2], 0, 2π - sep)]   # ∈ [0, 2π] × [0, 2π]
    end

    # select a distance metric
    d_t = PeriodicEuclidean([2π, 2π])  # distance function over a torus manifold

    # # offsets
    # offsets = [
    #     p -> torus_ψ1(p),
    #     p -> -torus_ψ1(p),
    #     p -> torus_ψ2(p),
    #     p -> -torus_ψ2(p),
    # ]

    # # one forms
    # Ω = OneForm[
    #     OneForm(1, (x, y) -> offset_size * torus_ψ1(x, y)),
    #     OneForm(1, (x, y) -> offset_size * -torus_ψ1(x, y)),
    #     OneForm(2, (x, y) -> offset_size * torus_ψ2(x, y)),
    #     OneForm(2, (x, y) -> offset_size * -torus_ψ2(x, y)),
    # ]

    return if cantype == :single 
        SingleCAN(
            "torus",
            cover,
            n,
            ξ_t,
            d_t,
            k_t;
            σ = :relu,
        )
    else
        CAN(
        "torus",
        cover,
        n,
        ξ_t,
        d_t,
        k_t;
        offset_size = offset_size,
        σ = :relu,
        α = 3.2,
        # offsets = offsets,
        # Ω = Ω
    )
    end
end


toruscan = torus_maker(:can)
toruscan_single = torus_maker(:single; k_t = LocalGlobalKernel(α = 1.0, σ = 50.0))