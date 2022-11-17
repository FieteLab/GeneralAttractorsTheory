using GeneralAttractors.Analysis
using Plots
import GeneralAttractors: load_data
using GeneralAttractors, Distances

"""
Create a PCA or Isomap visualization of a pointcloud
of neural activity from a simulation.
"""
ρ(x) = x
ρⁱ(x) = x

sim_fld = "abstract"
sim = "abstract_torus"
params = AnalysisParameters(
    max_nPC = 10,  # max num of PCs
    pca_pratio = 0.999999,       # fraction of variance explained
    n_isomap_dimensions = 3,
    isomap_k = 10,
    debug = true,   # avoid re-running analysis steps
)


# ------------------------- dimensionality reduction ------------------------- #
@info "PCA dimensionality reduction"
pca_dimensionality_reduction(sim_fld, sim, params)

@info "ISOMAP dimensionality reduction"
isomap_dimensionality_reduction(sim_fld, sim, params)


# ----------------------------------- plot ----------------------------------- #
X = load_data(sim_fld, sim*"_isomap_space")
scatter(eachrow(X)...)
