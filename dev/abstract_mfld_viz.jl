using GeneralAttractors.Analysis
using Plots
import GeneralAttractors: load_data
using GeneralAttractors, Distances
import GeneralAttractors.Analysis.ManifoldAnalysis: population_average
using Manifolds, Distances



ρ(x) = x
ρⁱ(x) = x

"""
Create a PCA or Isomap visualization of a pointcloud
of neural activity from a simulation.
"""



params = AnalysisParameters(
    max_nPC = 200,  # max num of PCs
    pca_pratio = 0.999999,       # fraction of variance explained
    n_isomap_dimensions = 3,
    isomap_k = 20,
    isomap_downsample = 10,
    debug = false,   # avoid re-running analysis steps
)

# folder where stuff is saved and the manifold/simulation name
sim_fld = "abstract"
mfld_name = "sphere"

# load data
simulations = []
for i in 1:50
    sim = "$(mfld_name)_$(i)_$(mfld_name)"
    history = load_simulation_history(sim_fld, sim*"_history")
    push!(simulations, population_average(history))
end
S = hcat(simulations...)


# ------------------------- dimensionality reduction ------------------------- #
# @info "PCA dimensionality reduction"
# pca_dimensionality_reduction(S, mfld_name, "mfld",  params; visualize=true)

# @info "ISOMAP dimensionality reduction"
# isomap_dimensionality_reduction(mfld_name, "mfld", params)


# ----------------------------------- plot ----------------------------------- #
X = load_data(mfld_name, "mfld_isomap_space")

# ? interactive plot
import GLMakie
fig = GLMakie.Figure(resolution=(1000,1000)); 
ax = GLMakie.Axis3(fig[1,1]); 
GLMakie.scatter!(ax,    eachrow(X)..., markersize=30, alpha=.5, strokecolor="black", color=:black)
fig


# ? static plot
scatter3d(
    eachrow(X)..., 
    markersize=10, alpha=.01, strokecolor="black", 
    camera=(45, 45), label=nothing, color=:black,
    xticks=nothing, xlabel="PC1",
    yticks=nothing, ylabel="PC2",
    zticks=nothing, zlabel="PC3",
    grid=false
    )

