using BSON
using ObjectivePaths
using NPZ
import MultivariateStats
import ManifoldLearning
import NearestNeighbors

DATADIR = path(pwd()) / "data"

""" 
    Check if a file at a path exists
"""
function checkpath(fld::String, savename::String, extension="bson")::Bool
    p = DATADIR / fld / "$savename.$extension"
    exists(p)
end


function savepath(fld::String, savename::String, extension="bson")
    # get destination path
    base = DATADIR / fld 
    !exists(base) && mkdir(base)
    svp = base / "$savename.$extension"
    exists(svp) && @warn "Saving $extension file at $(svp.path) - overwriting data"
    return svp.path
end

function loadpath(fld::String, loadname::String, extension="bson")
    # get destination path
    lp = DATADIR / fld / "$loadname.$extension"
    exists(lp) || error("Destination file not found at $(lp.path)")
    lp.path
end
# ---------------------------------------------------------------------------- #
#                                    SAVING                                    #
# ---------------------------------------------------------------------------- #
save_simulation_history(history, save_fld::String, savename::String) = bson(
        savepath(save_fld, savename), 
        Dict(:history=>history),
)
    
save_model(model, save_fld::String, savename::String, model_name::Symbol) = bson(
        savepath(save_fld, savename), 
        Dict(model_name=>model),
)


save_data(X::Array, save_fld::String, savename::String) =  npzwrite(
        savepath(save_fld, savename, "npz"), X
    )

# ---------------------------------------------------------------------------- #
#                                    LOADING                                   #
# ---------------------------------------------------------------------------- #
load_simulation_history(load_fld::String, loadname::String) = BSON.load(
    loadpath(load_fld, loadname)
)[:history]

load_model(load_fld::String, loadname::String, model_name::Symbol) = BSON.load(
    loadpath(load_fld, loadname)
)[model_name]

load_data(load_fld::String, loadname::String) = npzread(
    loadpath(load_fld, loadname, "npz")
)
