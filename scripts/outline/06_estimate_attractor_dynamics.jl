using NearestNeighbors

"""
Estimate attractor dynamics by:

1. loading torus simulation data
2. picking a sample plot
3. for n times, perturn the state in a random direction and measure
    dynamics along the manifold and off the manifold. 
4. to get the on/off manifold dynamics, use local PCA to get the tangent
    space and project the perturbation onto it.
"""

include("settings.jl")
move_to_datadir(supervisor, "mfld_top")


# load torus data
tag = "manifold_topology_data"
filters = Dict{Symbol, Any}(
    :tag => tag,
    :can => "torus",
)

@isdefined(X) || (X = load_and_concat_activations(; filters...) )
@isdefined(can) || (can = make_standard_single_torus_can())


# @time nntree = KDTree(X; reorder = false, leafsize = 5)