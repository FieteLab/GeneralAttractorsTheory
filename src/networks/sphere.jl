function sphere_maker(
    cantype;
    n::Int = 48,
    offset_size = 0.2,
    k = LocalGlobalKernel(α = 2.5, σ = 40.5),
    α = 46,
    σ = :softrelu, kwargs...
)

    # get neurons on S² ⊂ ℝ³
    X = fibonacci_sphere(n^2)

    # get neurons indices
    I = [(i,) for i = 1:size(X, 2)]

    # distance metric on the unit sphere
    d = SphericalDistance()
    
    # cover space
    cover = CoverSpace(S²)  # trivial cover space

    # define offset vector fields
    offsets = [
        p -> sphere_ψx(p), p -> -sphere_ψx(p),
        p -> sphere_ψy(p), p -> -sphere_ψy(p), 
        p -> sphere_ψz(p), p -> -sphere_ψz(p)
    ]

    # define one forms
    Ω = [
        OneForm(1, (p) -> offset_size * sphere_ψx(p)),
        OneForm(1, (p) -> -offset_size * sphere_ψx(p)),
        OneForm(2, (p) -> offset_size * sphere_ψy(p)),
        OneForm(2, (p) -> -offset_size * sphere_ψy(p)),
        OneForm(3, (p) -> offset_size * sphere_ψz(p)),
        OneForm(3, (p) -> -offset_size * sphere_ψz(p)),
    ]


    # construct CAN
    return if cantype == :single
        SingleCAN(
            "sphere",
            cover,
            (n, n),
            X,
            d,
            k;
            σ = σ,
        )
    else
        CAN(
            "sphere",
            cover,
            (n, n),
            I,
            X,
            d,
            k;
            offset_size = offset_size,
            offsets = offsets,
            Ω = Ω,
            α = α,
            σ = σ,
        )
    end
end

spherecan = sphere_maker(:can)
spherecan_single = sphere_maker(:single)