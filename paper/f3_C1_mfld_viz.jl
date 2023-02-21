"""
Isomap embedding visualization of activity manifold
"""

include("settings.jl")
move_to_datadir(supervisor, "mfld_top")
pyplot()

η = 1.0
_η = replace(string(η), "." => "_")
Δ = 1  # donwsampling

camera_angles = (100, 10),  (150, 70)

for (network, color) in zip(networks, networks_colors)
    (η > 0 && network ∉ ("ring", "torus", "sphere")) && continue

    # load data
    filters = Dict{Symbol, Any}(
        :tag => "d3_embeddings",
        :dim => 3,
        :can => network,
        :extension => "npz",
        :η => η
    )

    M = ProjectSupervisor.fetch(supervisor; filters...)[2][1]

    zlims = if network == "plane" 
        [-20, 20]
    elseif network == "mobius"
        [-10, 10]
    else
        [-20, 20]
    end

    plots = []
    for (i, camera) in enumerate(camera_angles)
        push!(plots, 
            scatter3d(
                M[1, 1:Δ:end], M[2, 1:Δ:end], M[3, 1:Δ:end],
                msa=0, msw=0,
                markersize = 4, # 1.0,
                alpha=.6,
                legend = false,
                title = network * " (η = $η)",
                xlabel = "Iso 1", ylabel = "Iso 2", zlabel = "Iso 3",
                color = color, 
                camera = camera,
                zlim = zlims,
                aspect_ratio=network == "line" ? :equal : :auto,
                showaxis = false,
                axis=nothing,
            )
        )
        
    end

    fig = plot(plots..., layout = (2, 1), size = (800, 800))
    save_plot(supervisor, fig, "f3_mfld_top_$(network)_noise_$(_η).png")
    # display(fig)
    # break
end

println("done")