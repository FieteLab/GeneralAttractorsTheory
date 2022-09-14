using GeneralAttractors
using GeneralAttractors.Analysis
using Plots

sim_name = "torus_sim"  # name of the simulation being anlyzed
# pca_dimensionality_reduction(sim_name)

# isomap_dimensionality_reduction(sim_name);

M = load_data(sim_name * "_iso")
# size(M)

animate_3d_scatter(M, sim_name*"_iso_anim"; alpha=.25, title="cacca")