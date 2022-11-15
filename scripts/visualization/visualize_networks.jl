using GeneralAttractors


# ----------------------------------- ring ----------------------------------- #
# show_connectivity(ring_attractor; plot_title = "Ring attractor connectivity")

# ----------------------------------- torus ---------------------------------- #
# include("../networks/torus.jl")
# show_connectivity(toruscan; plot_title = "Torus attractor connectivity")

# ---------------------------------- mobius ---------------------------------- #
include("../networks/mobius.jl")
show_connectivity(mobius_attractor; plot_title = "Mobius attractor connectivity")

