using GeneralAttractors.Analysis
using Plots
import GeneralAttractors: load_data
using GeneralAttractors, Distances
import GeneralAttractors.Analysis.ManifoldAnalysis: population_average, tda_on_pointcloud
using Manifolds, Distances, ManifoldsBase

# include("../scripts/networks/sphere.jl") 

ρ(x) = x
ρⁱ(x) = x

sim_fld = "sphere"
mfld_name = "decoding"


sim = "$(mfld_name)_sphere"
history = load_simulation_history(sim_fld, sim*"_history")

S = population_average(history)[:, 35]
th = .001
active = S .>= th
inactive = S .< th


sel = spherecan.X[2, :] .> 0.1

# scatter3d(
#     eachrow(
#     spherecan.X[:, inactive])...,
#     marker_z=S[inactive],
#     ms=3,
#     lw=2, 
#     msa=0, msw=0,
#     alpha=.2,
# )
scatter3d(
    eachrow(
    spherecan.X[:, sel .* inactive])...,
    marker_z=S[sel .* inactive],
    ms=3,
    lw=2, 
    # msa=0.0, 
    msw=0.1,
    # msc=S[sel],
    # alpha=(S./maximum(S) .+ .2)[sel],
    alpha=.3,
    xticks=nothing, yticks=nothing, zticks=nothing,
    camera=(180, 00),
    xlim=[-1, 1],
    ylim=[-1, 1],
    zlim=[-1, 1],
    clims=(0, .05),
    colorbar_title="activation",

)

scatter3d!(
    eachrow(
    spherecan.X[:, active])...,
    marker_z=S[active],
    ms=3,
    lw=2, 
    # msa=0.0, 
    msw=0.0,
    # msc=S[sel],
    xticks=nothing, yticks=nothing, zticks=nothing,
    camera=(180, 00),
    xlim=[-1, 1],
    ylim=[-1, 1],
    zlim=[-1, 1],
    clims=(0, .05),
    colorbar_title="activation",

)