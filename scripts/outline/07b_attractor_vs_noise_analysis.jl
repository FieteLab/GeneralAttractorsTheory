

"""
Load simulation data from the Torus with varying levelevs of noise and 
plot do TDA on them.
"""


include("settings.jl")
move_to_datadir(supervisor, "mfld_top_noise")


noise_values = (0.0, 0.5, 1.0)

SHOW_ISOMAP_3D = false
DO_TDA = true

# ------------------------------ plot 3d isomap ------------------------------ #
if SHOW_ISOMAP_3D
    plots = []
    for η in noise_values
        # load data
        filters = Dict{Symbol, Any}(
            :tag => "d3_embeddings",
            :can => "torus",
            :η => η,
        )

        meta, M = ProjectSupervisor.fetch(supervisor; filters...) 
        @assert length(M) == 1
        M = M[1]

        push!(plots, 
        scatter3d(
            M[1, 1:5:end], M[2, 1:5:end], M[3, 1:5:end],
            msa=0, msw=0,
            markersize = 3.0,
            legend = false,
            title = "η = $η",
            xlabel = "Iso 1", ylabel = "Iso 2", zlabel = "Iso 3",
            color = :black, 
            camera = (70, 25),
            alpha=.7,
            zlim = [-20, 20],
            showaxis = false,
            axis=nothing,
        )
        )
    end
    image.png
    fig = plot(plots..., 
            layout=(1, 3), 
            size=(1200, 800);
            plot_font_size_kwargs...
            )
    display(fig)
    save_plot(supervisor, fig, "07b_noise_isomap_3d")
end


# ------------------------------------ tda ----------------------------------- #

if DO_TDA
    for η in noise_values
        # load data
        filters = Dict{Symbol, Any}(
            :tag => "d10_embeddings",
            :can => "torus",
            :η => η,
        )

        do_tda(supervisor, filters, "07b_TDA_η=$η"; max_d=2)
    end
end


