using BSON
using ObjectivePaths

function save_simulation_history(history, savename::String)
    # get destination path
    savefld = path(pwd()) / "data"
    savepath = savefld / "$savename.bson"
    exists(savepath) && @warn "Saving history at $(savepath.path) - overwriting data"

    
    # save
    bson(
        savepath.path, 
        Dict(:history=>history),
    )
    println(savepath)

end


function load_simulation_history(savename::String)
    # get destination path
    savefld = path(pwd()) / "data"
    savepath = savefld / "$savename.bson"
    exists(savepath) || error("Destination file not found at $(savepath.path)")

    BSON.load(savepath.path)[:history]
end
