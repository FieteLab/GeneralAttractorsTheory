
function mobius_maker(
    cantype; 
    offset_size = 0.2,
    n = nothing,  # not used but here for consistency with other methods
    α = 35,
    σ = :softrelu,
    k = LocalGlobalKernel(α = 2.5, σ = 1.5)
)
    mfld = Mobius()

    # number of neurons
    y_extent = abs(mfld.xmin[1]) + abs(mfld.xmax[1])
    n = isnothing(n) ? 
            ((Int ∘ round)(y_extent / 0.025), (Int ∘ round)(2π / 0.1)) :
            (n, n)

    # cover space
    cover = CoverSpace(mfld)

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
        OneForm(1, (t, θ) -> offset_size * MB_ψ1(t, θ)),
        OneForm(2, (t, θ) -> -offset_size * MB_ψ1(t, θ)),
        OneForm(3, (t, θ) -> offset_size * MB_ψ2(t, θ)),
        OneForm(4, (t, θ) -> -offset_size * MB_ψ2(t, θ)),
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