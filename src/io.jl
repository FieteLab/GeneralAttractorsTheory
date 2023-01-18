import DrWatson: datadir, savename, tagsave, wsave

function save_sim_data(folder_name::String, savename_data::Dict,  sim, history, X)
    save_path = datadir("simulations", folder_name, savename("metadata", savename_data, "bson"))
    tagsave(save_path, savename_data)

    save_path = datadir("simulations", folder_name, savename("history", savename_data, "bson"))
    wsave(save_path, Dict(:history => history))

    save_path = datadir("simulations", folder_name, savename("trajectory", savename_data, "npz"))
    wsave(save_path, sim.trajectory.X)

    save_path = datadir("simulations", folder_name, savename("decoded_X", savename_data, "npz"))
    wsave(save_path, X)
end