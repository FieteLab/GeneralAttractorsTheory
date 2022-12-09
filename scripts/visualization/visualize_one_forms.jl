using Plots
using Distances

using GeneralAttractors
using GeneralAttractors.Kernels
using GeneralAttractors: lerp
using GeneralAttractors.ManifoldUtils

include("../networks/torus.jl")

# ----------------------------------- show ----------------------------------- #
dx, scale = 8, 2

i = 3
p = show_oneforms(
    toruscan.Ω[3],
    toruscan.C,
    cover.M.xmin,
    cover.M.xmax;
    dx = dx,
    scale = scale,
    color = :red,
)

# plots = []
# for i = 1:2:length(toruscan.Ω)
#     p = show_oneforms(
#         toruscan.Ω[i],
#         toruscan.C,
#         cover.M.xmin,
#         cover.M.xmax;
#         dx = dx,
#         scale = scale,
#         color = :red,
#     ) |> display
#     break
#     # push!(plots, p)
# end
# # plot(plots...)
