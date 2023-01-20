using Term.Progress
include("settings.jl")

"""
Load data with single CAN simulations for each network and plot the 3d embedding.
"""

move_to_datadir(supervisor, "mfld_top")


plots = []
camera_angles = (45, 25), (20, 20), (70, 70)

for network in networks
    # load data
    filters = Dict{Symbol, Any}(
        :tag => "d3_embeddings",
        :can => network,
        :dim => 3,
    )

    _, M = ProjectSupervisor.fetch(supervisor; filters...) 
    @assert length(M) == 1
    M = M[1]
    @assert size(M, 1) == 3


    for camera in camera_angles
        push!(plots, 
            scatter3d(
                M[1, 1:25:end], M[2, 1:25:end], M[3, 1:25:end],
                msa=0, msw=0,
                markersize = 4.0,
                legend = false,
                title = network,
                xlabel = "Iso 1", ylabel = "Iso 2", zlabel = "Iso 3",
                color = :black, 
                camera = camera,
                alpha=.2,
                # showaxis = false,
                # axis=nothing,
            )
        )
        break
    end

end

fig = plot(plots..., 
        layout=(3, 2), 
        size=(800, 800)
        )

