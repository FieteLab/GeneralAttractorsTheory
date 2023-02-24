function sphere_maker(
    cantype;
    n::Int = 48,
    offset_size = 0.2,
    use_offset_fields = false,
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
    normalize(x) = norm(x) > 0 ? x ./ norm(x) : x
    ψ_x = use_offset_fields ? normalize ∘ sphere_ψx : sphere_ψx
    ψ_y = use_offset_fields ? normalize ∘ sphere_ψy : sphere_ψy
    ψ_z = use_offset_fields ? normalize ∘ sphere_ψz : sphere_ψz

    offsets = [
        p -> ψ_x(p), p -> -ψ_x(p),
        p -> ψ_y(p), p -> -ψ_y(p), 
        p -> ψ_z(p), p -> -ψ_z(p)
    ]

    # define one forms
    Ω = [
        OneForm(1, (p) -> offset_size * ψ_x(p)),
        OneForm(1, (p) -> -offset_size * ψ_x(p)),
        OneForm(2, (p) -> offset_size * ψ_y(p)),
        OneForm(2, (p) -> -offset_size * ψ_y(p)),
        OneForm(3, (p) -> offset_size * ψ_z(p)),
        OneForm(3, (p) -> -offset_size * ψ_z(p)),
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

# spherecan = sphere_maker(:can)
# spherecan_single = sphere_maker(:single)