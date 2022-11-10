using GeneralAttractors

include("../networks/torus.jl")

# show_connectivity(ring_attractor; plot_title = "Ring attractor connectivity")

show_connectivity(toruscan; plot_title = "Torus attractor connectivity")

# show_connectivity(mobius_attractor; plot_title = "Mobius attractor connectivity")

# show_connectivity(
#     sph;
#     xlabel = "longitude",
#     ylabel = "latitude",
#     plot_title = "Sphere attractor connectivity",
#     idxs = [1, 30, 2100, 1500, 800],
#     aspect_ratio=0.5,
#     size=(800, 400)
# )
