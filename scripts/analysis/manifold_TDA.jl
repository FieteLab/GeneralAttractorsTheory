using GeneralAttractors.Analysis
using Plots
import GeneralAttractors: load_data
using GeneralAttractors, Distances
import GeneralAttractors.Analysis.ManifoldAnalysis: population_average, tda_on_pointcloud
using Manifolds, Distances, ManifoldsBase

ρ(x) = x
ρⁱ(x) = x

params = AnalysisParameters(
    max_nPC = 10,  # max num of PCs
    pca_pratio = 0.999999,       # fraction of variance explained
    n_isomap_dimensions = 3,
    isomap_k = 20,
    isomap_downsample = 20,
    tda_threshold = 10,       # threshold to reduce TDA computation 
    tda_downsample_factor = 100,        # temporal downsampling of data for TDA
    tda_dim_max = 2,        # max feature dimension, starting at 0
    debug = true,   # avoid re-running analysis steps
)

# folder where stuff is saved and the manifold/simulation name
sim_fld = "abstract"
mfld_name = "torus"
n_sims = 50

# load data
simulations = []
for i = 1:n_sims
    sim = "$(mfld_name)_$(i)_$(mfld_name)"
    history = load_simulation_history(sim_fld, sim * "_history")
    push!(simulations, population_average(history))
end
S = hcat(simulations...)




# ------------------------- dimensionality reduction ------------------------- #
@info "PCA dimensionality reduction"
pca_dimensionality_reduction(S, mfld_name, "mfld", params; visualize = true)


# ------------------------------------ TDA ----------------------------------- #

tda_on_pointcloud(mfld_name, "mfld", params)
