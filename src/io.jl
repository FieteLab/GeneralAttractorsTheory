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



