using Plots
using Distances

using GeneralAttractors
using GeneralAttractors.Kernels
using GeneralAttractors: lerp
using GeneralAttractors.Manifolds



cover = CoverSpace(ℝ², T, (x, y) -> [mod(x - 32, 64), mod(y - 32, 64)])



n = (64, 64)
function ξ_t(i::Int, j::Int)::Vector  # neurons coordinates function
    n̂_i, n̂_j = Int(n[1] / 2), Int(n[2] / 2)
    [lerp(i, n[1], -n̂_i, n̂_i), lerp(j, n[2], -n̂_j, n̂_j)]   # ∈ [-n/2, n/2] × [-n/2, n/2]
end
d_t = PeriodicEuclidean([n...])  # distance function over a torus manifold

# connectivity kernel 
k_t = DiffOfExpKernel(; λ = 13.0)



Ω = OneForm[
    OneForm(1, x -> sin(x / n[1]) + 1.25),
    OneForm(1, x -> -(sin(x / n[1]) + 1.25)),
    OneForm(2, x -> sin(x / n[2]) + 1.25),
    OneForm(2, x -> -(sin(x / n[2]) + 1.25)),
]


tor = CAN("torus", cover, n, ξ_t, d_t, k_t; Ω = Ω)


p = show_oneforms(tor.Ω[1], tor.C, [-100, -100], [100, 100]; dx = 10, scale = 1)
show_oneforms!(
    p,
    tor.Ω[2],
    tor.C,
    [-100, -100],
    [100, 100];
    dx = 10,
    scale = 1,
    color = :red,
)
show_oneforms!(
    p,
    tor.Ω[3],
    tor.C,
    [-100, -100],
    [100, 100];
    dx = 10,
    scale = 1,
    color = :green,
)
show_oneforms!(
    p,
    tor.Ω[4],
    tor.C,
    [-100, -100],
    [100, 100];
    dx = 10,
    scale = 1,
    color = :blue,
)
plot!(
    [-32, 32],
    [-32, -32],
    lw = 4,
    ls = :dash,
    color = :black,
    label = nothing,
    alpha = 0.25,
)
plot!(
    [-32, 32],
    [32, 32],
    lw = 4,
    ls = :dash,
    color = :black,
    label = nothing,
    alpha = 0.25,
)
plot!(
    [-32, -32],
    [-32, 32],
    lw = 4,
    ls = :dash,
    color = :black,
    label = nothing,
    alpha = 0.25,
)
plot!(
    [32, 32],
    [-32, 32],
    lw = 4,
    ls = :dash,
    color = :black,
    label = nothing,
    alpha = 0.25,
)


# TODO get well defined behavior for ρ such that one-forms act correctly
# TODO test: implement decoding
