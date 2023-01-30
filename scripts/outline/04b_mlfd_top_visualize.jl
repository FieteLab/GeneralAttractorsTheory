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
    )

    meta, M = ProjectSupervisor.fetch(supervisor; filters...) 
    @assert length(M) == 2 network
    M = M[2]
    @assert size(M, 1) == 3 network


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
                markersize = 3.0,
                legend = false,
                title = network,
                xlabel = "Iso 1", ylabel = "Iso 2", zlabel = "Iso 3",
                color = color, 
                camera = camera,
                alpha=.5,
                zlim = zlims,
                showaxis = false,
                axis=nothing,
            )
        )
        # break
    end

end

fig = plot(plots..., 
        layout=(length(networks), 2), 
        size=(800, 1000)
        )

save_plot(supervisor, fig, "04b_mfld_top_pointclouds")
