using Plots
using Distances

using GeneralAttractors
using GeneralAttractors.Kernels
using GeneralAttractors: lerp
using GeneralAttractors.Manifolds

# ------------------------------- make network ------------------------------- #

n = (64, 64)
function ξ_s(i::Int, j::Int)::Vector
    [lerp(i, n[1], -π, π), lerp(j, n[2], -π / 2, π / 2)]
end
d_s = SphericalAngle()
k_s = DiffOfExpKernel(; λ = 0.75)

cover = CoverSpace(S², S², (x, y) -> [x, y])


Ω = OneForm[
    OneForm(1, x -> 1),
    OneForm(1, x -> -1),
    OneForm(2, x -> 1),
    OneForm(2, x -> -1),
]

can = CAN("sphere", cover, n, ξ_s, d_s, k_s; offset_size=[0.2, 0.2, 0.1, 0.1], Ω=Ω)


# ----------------------------------- show ----------------------------------- #
dx, scale = 1, 0.2
p = show_oneforms(can.Ω[1], can.C, cover.M.xmin, cover.M.xmax; dx = dx, scale = scale)
show_oneforms!(
    p,
    can.Ω[2],
    can.C,
    cover.M.xmin,
    cover.M.xmax;
    dx = dx,
    scale = scale,
    color = :red,
)
show_oneforms!(
    p,
    can.Ω[3],
    can.C,
    cover.M.xmin,
    cover.M.xmax;
    dx = dx,
    scale = scale,
    color = :green,
)
show_oneforms!(
    p,
    can.Ω[4],
    can.C,
    cover.M.xmin,
    cover.M.xmax;
    dx = dx,
    scale = scale,
    color = :blue,
)
