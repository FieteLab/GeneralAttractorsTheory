include("settings.jl")

"""
Load data with single CAN simulations for each network and plot the 3d embedding.
"""

move_to_datadir(supervisor, "mfld_top")


plots = []
camera_angles = (100, 25),  (70, 60)

for (network, color) in zip(networks, networks_colors)
    # load data
    filters = Dict{Symbol, Any}(
        :tag => "d3_embeddings",
        :can => network,
        :extension => "npz"
    )

    M = ProjectSupervisor.fetch(supervisor; filters...)[2][1]
    zlims = if network == "plane" 
        [-20, 20]
    elseif network == "mobius"
        [-10, 10]
    else
        [-20, 20]
    end

    for (i, camera) in enumerate(camera_angles)
        # i == 1 && continue
        push!(plots, 
            scatter3d(
                M[1, 1:5:end], M[2, 1:5:end], M[3, 1:5:end],
                msa=0, msw=0,
                markersize = 1.0,
                alpha=.25,
                legend = false,
                title = network,
                xlabel = "Iso 1", ylabel = "Iso 2", zlabel = "Iso 3",
                color = color, 
                camera = camera,
                zlim = zlims,
            #     showaxis = false,
            #     axis=nothing,
            )
        )
        # break
    end

end

fig = plot(plots..., 
        layout=(length(networks), 2), 
        size=(1200, 1400)
        )

save_plot(supervisor, fig, "04b_mfld_top_pointclouds")
fig     