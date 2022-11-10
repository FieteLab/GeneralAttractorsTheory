using GeneralAttractors


# ----------------------------------- ring ----------------------------------- #
# show_connectivity(ring_attractor; plot_title = "Ring attractor connectivity")

# ----------------------------------- torus ---------------------------------- #
# include("../networks/torus.jl")
# show_connectivity(toruscan; plot_title = "Torus attractor connectivity")

# ---------------------------------- mobius ---------------------------------- #
# show_connectivity(mobius_attractor; plot_title = "Mobius attractor connectivity")

# ---------------------------------- sphere ---------------------------------- #
include("../networks/sphere.jl")
show_connectivity(
    spherecan;
    xlabel = "longitude",
    ylabel = "latitude",
    plot_title = "Sphere attractor connectivity",
    # idxs = [1, 30, 2100, 1500, 800],
    aspect_ratio=1,
    size=(1000, 1000)
)
