using BSON
using ObjectivePaths
using NPZ
import MultivariateStats
import ManifoldLearning
import NearestNeighbors

DATADIR = path(pwd()) / "data"

function savepath(savename::String, extension="bson")
    # get destination path
    svp = DATADIR / "$savename.$extension"
    exists(svp) && @warn "Saving $extension file at $(svp.path) - overwriting data"
    return svp.path
end

function loadpath(loadname::String, extension="bson")
    # get destination path
    lp = DATADIR / "$loadname.$extension"
    exists(lp) || error("Destination file not found at $(lp.path)")
    lp.path
end
# ---------------------------------------------------------------------------- #
#                                    SAVING                                    #
# ---------------------------------------------------------------------------- #
save_simulation_history(history, savename::String) = bson(
        savepath(savename), 
        Dict(:history=>history),
)
    
save_model(model, savename::String, model_name::Symbol) = bson(
        savepath(savename), 
        Dict(model_name=>model),
)


save_data(X::Array, savename::String) =  npzwrite(
        savepath(savename, "npz"), X
    )

# ---------------------------------------------------------------------------- #
#                                    LOADING                                   #
# ---------------------------------------------------------------------------- #
load_simulation_history(loadname::String) = BSON.load(
    loadpath(loadname)
)[:history]

load_model(loadname::String, model_name::Symbol) = BSON.load(
    loadpath(loadname)
)[model_name]

load_data(loadname::String) = npzread(
    loadpath(loadname, "npz")
)
