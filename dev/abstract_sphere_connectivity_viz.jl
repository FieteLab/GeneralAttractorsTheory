using GeneralAttractors.Analysis
using Plots
import GeneralAttractors: load_data
using GeneralAttractors, Distances
import GeneralAttractors.Analysis.ManifoldAnalysis: population_average, tda_on_pointcloud
using Manifolds, Distances, ManifoldsBase

# include("../scripts/networks/sphere.jl") 


idx = 4000
W = spherecan.Ws[1][idx, :]


sel = spherecan.X[2, :] .< 0


scatter3d(
    eachrow(
    spherecan.X[:, sel])...,
    marker_z=W[sel],
    ms=3,
    lw=2, 
    # msa=0.0, 
    msw=0.1,
    # msc=S[sel],
    xticks=nothing, yticks=nothing, zticks=nothing,
    camera=(180, 00),
    xlim=[-1, 1],
    ylim=[-1, 1],
    zlim=[-1, 1],
    # clims=(0, .05),
    colorbar_title="weight",
    color=:bwr,

)

scatter3d!(
    [spherecan.X[1, idx]],
    [spherecan.X[2, idx]],
    [spherecan.X[3, idx]],
    ms=10, color=:black,
)