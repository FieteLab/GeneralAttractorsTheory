function torus_maker(cantype; 
        n::Int=48,
        k::AbstractKernel = LocalGlobalKernel(α = 2.5, σ = 5.0), 
        offset_size::Number = 0.2,
        α = 3.2,
        σ = :softrelu,
        use_offset_fields::Bool = false,
        cover_manifold::Symbol = :default,
    )
    # number of neurons
    n = (n, n) # number of neurons per dimension


    ρ(x, y) = if cover_manifold == :default
            [mod(x, 2π), mod(y, 2π)] 
    elseif cover_manifold == :cylinder
        [x, mod(y, 2π)] 
    end
    ρ(v) = ρ(v...)

    """
        ρⁱ(x, y)

    Inverse of the cover map rho over the domain.
    Given a point (x,y) in N it gives a set of points (x̂, ŷ)
    in the cover space such that ρ(x̂, ŷ)=(x,y)
    """
    function ρⁱ(x, y; n = 40)
        if cover_manifold == :default
            pts = zeros(2, n^2)
            for (c, i) in enumerate(-n/2:(n/2-1)), (k, j) in enumerate(-n/2:(n/2-1))
                x̂ = x + 2π * i
                ŷ = y + 2π * j
                pts[:, (c-1)*n+k] = [x̂, ŷ]
            end
        else
            pts = zeros(2, n^2)
            for (k, j) in enumerate(-n/2:(n/2-1))
                ŷ = y + 2π * j
                pts[:, k] = [x, ŷ]
            end
        end
        return pts
    end
    ρⁱ(w::Vector; n = 6) = ρⁱ(w...; n = n)
    
    cover = if cover_manifold == :default
        CoverSpace(Manifoldℝ²(25), Torus(), ρ, ρⁱ)
    elseif cover_manifold == :cylinder
        CoverSpace(Cylinder(20), Torus(), ρ, ρⁱ)
    end

    # define a function to get the coordinates of each neuron in the lattice
    function ξ(i::Int, j::Int)::Vector  # neurons coordinates function
        sep = 2π / n[1]
        [lerp(i, n[1], 0, 2π - sep), lerp(j, n[2], 0, 2π - sep)]   # ∈ [0, 2π] × [0, 2π]
    end

    # select a distance metric
    d = PeriodicEuclidean([2π, 2π])  # distance function over a torus manifold

    # offsets
    if use_offset_fields
        offsets = [
            p -> torus_ψ1(p),
            p -> -torus_ψ1(p),
            p -> torus_ψ2(p),
            p -> -torus_ψ2(p),
        ]

        # one forms
        Ω = OneForm[
            OneForm(1, (v) -> offset_size * torus_ψ1(v)),
            OneForm(1, (v) -> offset_size * -torus_ψ1(v)),
            OneForm(2, (v) -> offset_size * torus_ψ2(v)),
            OneForm(2, (v) -> offset_size * -torus_ψ2(v)),
        ]
    else
        offsets, Ω = nothing, nothing
    end

    return if cantype == :single 
        SingleCAN(
            "torus",
            cover,
            n,
            ξ,
            d,
            k;
            σ = σ,
        )
    else
        CAN(
        "torus",
        cover,
        n,
        ξ,
        d,
        k;
        offset_size = offset_size,
        σ = σ,
        α = α,
        offsets = offsets,
        Ω = Ω
    )
    end
end


# toruscan = torus_maker(:can)
# toruscan_single = torus_maker(:single; k = LocalGlobalKernel(α = 1.0, σ = 50.0))