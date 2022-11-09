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

O = [
    [cos(0), sin(0)],
    [cos(π * 1/4), sin(π * 1/4)],
    [cos(π * 2/4), sin(π * 2/4)],
    [cos(π * 3/4), sin(π * 3/4)],
    [cos(π), sin(π)],
    [cos(π+π*1/4), sin(π+π*1/4)],
    [cos(π+π*2/4), sin(π+π*2/4)],
    [cos(π+π*3/4), sin(π+π*3/4)],
]

Ω = OneForm[
    OneForm(1, (x, y) -> O[1]),
    OneForm(1, (x, y) -> O[2]),
    OneForm(1, (x, y) -> O[3]),
    OneForm(1, (x, y) -> O[4]),
    OneForm(2, (x, y) -> O[5]),
    OneForm(2, (x, y) -> O[6]),
    OneForm(2, (x, y) -> O[7]),
    OneForm(2, (x, y) -> O[8]),
]


can = CAN("sphere", cover, n, ξ_s, d_s, k_s; 
        offset_size=0.15,
        φ=sphere_embedding,
        offsets=O,
        Ω=Ω
        )

# ----------------------------------- show ----------------------------------- #
dx, scale = 1, 0.2

plt = plot()
for i in 1:length(Ω)
    show_oneforms!(
        plt,
        can.Ω[i],
        can.C,
        cover.M.xmin,
        cover.M.xmax;
        dx = dx,
        scale = scale,
        color = :red,
    )
end
plt