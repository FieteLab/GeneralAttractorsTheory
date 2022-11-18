using GeneralAttractors.Analysis
using Plots
import GeneralAttractors: load_data
using GeneralAttractors, Distances
import GeneralAttractors.Analysis.ManifoldAnalysis: population_average, tda_on_pointcloud
using Manifolds, Distances, ManifoldsBase

ρ(x) = x
ρⁱ(x) = x

params = AnalysisParameters(
    max_nPC = 400,  # max num of PCs
    pca_pratio = 0.999999,       # fraction of variance explained
    tda_threshold = 0.5,       # threshold to reduce TDA computation 
    tda_downsample_factor = 20,        # temporal downsampling of data for TDA
    tda_dim_max = 1,        # max feature dimension, starting at 0

    debug = true,   # avoid re-running analysis steps
)

# folder where stuff is saved and the manifold/simulation name
sim_fld = "abstract"
mfld_name = "sphere"
n_sims = 50

# load data
simulations = []
for i in 1:n_sims
    sim = "$(mfld_name)_$(i)_$(mfld_name)"
    history = load_simulation_history(sim_fld, sim*"_history")

    s = population_average(history)
    # i = round(size(s, 2), digits=-1) |> Int
    # s = s[:, 1:i]
    # n = size(s, 1)
    # s = mean(reshape(s, n, 5, :), dims=2)
    # s = reshape(s, n, :)
    
    push!(simulations, s)
end
S = hcat(simulations...)



# ------------------------- dimensionality reduction ------------------------- #
@info "PCA dimensionality reduction"
pca_dimensionality_reduction(S, mfld_name, "mfld",  params; visualize=false)


# ------------------------------------ TDA ----------------------------------- #

tda, plt = tda_on_pointcloud(mfld_name, "mfld", params)
plt