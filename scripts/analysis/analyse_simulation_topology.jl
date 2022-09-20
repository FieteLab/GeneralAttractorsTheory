using Plots
using Term
using Statistics

using GeneralAttractors
using GeneralAttractors.Analysis

"""
Load data from a CAN simulation and analyze
the topology of the activity manifold.
"""

sim = "torus_sim"
params = AnalysisParameters(
    intrinsic_d_neighborhood_size=50,
    debug=true,   # avoid re-running analysis steps
)
@info "Running manifold analysis for simulation '$sim'"
tprintln(params)

# ------------------------- dimensionality reduction ------------------------- #
# pca_dimensionality_reduction(sim, params)
isomap_dimensionality_reduction(sim, params)


# ------------------------- intrinsic dimensionality ------------------------- #
d = estimate_intrinsic_dimensionality(sim, params)
μ, σ = round(mean(d); digits=3), round(std(d); digits=2)
print(Panel("Intrinsic dimensionality: $μ ± $σ", title="Local PCA", style="green"))

# ------------------------------------ TDA ----------------------------------- #
# estimate_manifold_topology(sim, params)