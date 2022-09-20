using Plots
using Term

using GeneralAttractors
using GeneralAttractors.Analysis

"""
Load data from a CAN simulation and analyze
the topology of the activity manifold.
"""

sim = "torus_sim"
params = AnalysisParameters(
    debug=false,   # avoid re-running analysis steps
)
@info "Running manifold analysis for simulation '$sim'"
println(params)

# ------------------------- dimensionality reduction ------------------------- #
pca_dimensionality_reduction(sim, params)
isomap_dimensionality_reduction(sim, params)


# ------------------------- intrinsic dimensionality ------------------------- #
d = estimate_intrinsic_dimensionality(sim, params)
print(Panel("Intrinsic dimensionality: $d", style="green"))

# ------------------------------------ TDA ----------------------------------- #
estimate_manifold_topology(sim, params)