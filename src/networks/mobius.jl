
function mobius_maker(
    cantype; 
    offset_size = 0.2,
    n = nothing,  # not used but here for consistency with other methods
    α = 35,
    σ = :softrelu,
    k = LocalGlobalKernel(α = 2.5, σ = 1.5), 
    cover_manifold::Symbol = :default,
    kwargs...
)
    mfld = Mobius()

    # number of neurons
    y_extent = abs(mfld.xmin[1]) + abs(mfld.xmax[1])
    n = isnothing(n) ? 
            ((Int ∘ round)(y_extent / 0.02), (Int ∘ round)(2π / 0.2)) :
            (n, n)

    # cover space
    if cover_manifold == :default
        cover = CoverSpace(mfld)
    else
        # create a cover space from the cylinder to the mobius strip
        ρ(x, y) = [y, x]
        ρ(v) = ρ(v...)

        ρⁱ(x, y) = [y, x]
        ρⁱ(v) = ρⁱ(v...)


        cover = CoverSpace(
            Cylinder(y_extent/2), mfld, ρ, ρⁱ;
        )
    end

    # coordinates function (from neurons index to lattice coordintes)
    sep = 2π / n[2]
    ξ(i::Int, j::Int)::Vector = [
        lerp(i, n[1], mfld.xmin[1], mfld.xmax[1]),
        lerp(j, n[2], 0, 2π - sep),
    ]

    # metric
    d = MobiusEuclidean()

    # define offset vector fields
    offsets = [
        p -> MB_ψ1(p),
        p -> -MB_ψ1(p),
        p -> MB_ψ2(p),
        p -> -MB_ψ2(p),
    ]

    # define one forms
    Ω = [
        OneForm(1, (p) -> offset_size * MB_ψ1(p)),
        OneForm(1, (p) -> -offset_size * MB_ψ1(p)),
        OneForm(2, (p) -> offset_size * MB_ψ2(p)),
        OneForm(2, (p) -> -offset_size * MB_ψ2(p)),
    ]

    # construct network
    return if cantype == :single
        SingleCAN(
            "mobius",
            cover,
            n,
            ξ,
            d,
            k;
            σ = σ,
        )
    else
        CAN(
            "mobius",
            cover,
            n,
            ξ,
            d,
            k;
            offset_size = offset_size,
            offsets = offsets,
            Ω = Ω,
            σ = σ,
            α = α,
        )

    end

end


mobiuscan = mobius_maker(:can)
mobiuscan_single = mobius_maker(:single)