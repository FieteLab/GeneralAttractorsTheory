using GeneralAttractors.Analysis
using Plots
import GeneralAttractors: load_data
using GeneralAttractors, Distances
import GeneralAttractors.Analysis.ManifoldAnalysis: population_average


ρ(x) = x
ρⁱ(x) = x

"""
Create a PCA or Isomap visualization of a pointcloud
of neural activity from a simulation.
"""



params = AnalysisParameters(
    max_nPC = 500,  # max num of PCs
    pca_pratio = 0.999999,       # fraction of variance explained
    n_isomap_dimensions = 3,
    isomap_k = 20,
    debug = false,   # avoid re-running analysis steps
)


sim_fld = "abstract"


simulations = []
for i in 1:20
    sim = "torus_$(i)_torus"
    history = load_simulation_history(sim_fld, sim*"_history")
    push!(simulations, population_average(history))
end
S = hcat(simulations...)


# # ------------------------- dimensionality reduction ------------------------- #
@info "PCA dimensionality reduction"
# pca_dimensionality_reduction(S,sim_fld, "torus",  params; visualize=true)

@info "ISOMAP dimensionality reduction"
# isomap_dimensionality_reduction(sim_fld, "torus", params)


# # ----------------------------------- plot ----------------------------------- #
X = load_data(sim_fld, "torus_isomap_space")
# X = load_data(sim_fld, sim*"_pca_space")

import GLMakie
GLMakie.scatter(eachrow(X)...)
