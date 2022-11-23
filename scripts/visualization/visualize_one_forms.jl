using Plots
using Distances

using GeneralAttractors
using GeneralAttractors.Kernels
using GeneralAttractors: lerp
using GeneralAttractors.ManifoldUtils

include("../networks/sphere.jl")

# ----------------------------------- show ----------------------------------- #
dx, scale = 0.5, 0.25

plots = []
for i = 1:2:length(Ω)
    p = show_oneforms(
        spherecan.Ω[i],
        spherecan.C,
        cover.M.xmin,
        cover.M.xmax;
        dx = dx,
        scale = scale,
        color = :red,
    )
    push!(plots, p)
end
plot(plots...)
