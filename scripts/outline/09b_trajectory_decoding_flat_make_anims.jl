using Plots

using Term.Progress
using GeneralAttractors
using GeneralAttractors.Simulations
import GeneralAttractors.Simulations: plot_trajectory_and_decoded
import GeneralAttractors: animate_simulation_data

using MultivariateStats, ManifoldLearning

include("settings.jl")
tag = "decoding_data"


# ---------------------------- plotting functions ---------------------------- #

spatial_downsampling = 5
temporal_downsampling = 5
fps = 25

for network in ("plane", "cylinder", "torus")
    # make network to get neurons lattice coordinates
    can = network_makers[network](:single; n=48)
    w_x = range(can.X[1,1], can.X[1,end]; length=can.n[1])
    w_y = range(can.X[2,1], can.X[2,end]; length=can.n[2])

    # --------------------------------- get data --------------------------------- #
    # load 3d iso embedding model and embedded data
    set_datadir(supervisor,  (supervisor.datadir / "mfld_top").path)
    embeddings, M = ProjectSupervisor.fetch(
        supervisor; 
        tag="d3_embeddings", 
        can=network, 
        dim=3
    )[2]
    pca = embeddings[:pca]
    iso = embeddings[:iso]
    @assert size(M, 1) == 3 size(M)


    # load path integration data
    set_datadir(supervisor,  ((supervisor.datadir-1) / "path_int").path)
    filters = Dict{Symbol, Any}(
        :tag => tag,
        :network => network,
        :funky => false,
        :duration => 2500,
    )

    data = ProjectSupervisor.fetch(supervisor; filters...)[2][1]
    history = data["h"]
    trajectory = data["trajectory"]
    S = sum(history.S; dims=2)[:, 1, :]
    S_embedd = predict(iso, predict(pca, S))
    X = history.x_M_decoded


    @info "Got data" S M iso

    # ------------------------------ make animation ------------------------------ #

    nframes = size(trajectory.X, 1)
    skipframes = trajectory.still / dt |> floor |> Int

    # compute fps
    n_anim_frames = (nframes - skipframes)/temporal_downsampling  |> round |> Int
    @info "Ready to animate" nframes skipframes n_anim_frames fps spatial_downsampling temporal_downsampling

    anim = Animation()
    pbar = ProgressBar()
    Progress.with(pbar) do
        job = addjob!(pbar, description = "Making animation: $network", N = n_anim_frames)
        for i in 1:temporal_downsampling:nframes
            i < skipframes && continue
            fnum = max(i - skipframes, 1)
            hist_frame = (fnum)  ÷ history.average_over + 1

            # plot trajectory & decoded
            p1 = plot(trajectory, fnum+skipframes)
            plot!(
                p1, X[1, 1:fnum], X[2, 1:fnum],
                lw = 6, color=:red, alpha=.5,
                label = "decoded",
                xlabel = "m₁", ylabel = "m₂",
                )

            # plot neurons activation on lattice
            s =  S[1:spatial_downsampling:end, hist_frame]
            p2 = scatter3d(
                can.X[1, 1:spatial_downsampling:end], 
                can.X[2, 1:spatial_downsampling:end], 
                s,
                marker_z = s, msa=0, msw=0,
                colorbar = false,
                xlabel = "θ₁", ylabel = "θ₂", zlabel = "",
                zlim=[0, 0.75], grid=false,
                clims=(-0.2, 0.75),
                label=nothing,
            )

            # plot neurons activation in ISO space
            p3 = scatter(
                eachrow(M[:, 1:10:end])...,
                color=:black, alpha=.5, ms=1, 
                showaxis = false,
                axis=nothing,
                xlabel = "ISO 1", ylabel = "ISO 2", zlabel = "ISO 3",
                label=nothing,
            )

            s_embed = S_embedd[
                :, 
                max(hist_frame-100, 1):hist_frame,
            ]
            scatter!(
                eachrow(s_embed)...,
                marker_z = 1:size(s_embed, 2),
                msa=0, msw=0, colorbar = false,
                label=nothing
            )

            plt = plot(p1, p2, p3, layout=(1,3), size=(1200, 800),
                right_margin = 12Plots.mm,
                left_margin = 12Plots.mm,
                top_margin = 12Plots.mm,
                bottom_margin = 12Plots.mm,
            )
            frame(anim)
            update!(job)
        end
    end

    gif(anim, (supervisor.projectdir / "plots" / "path_int_$(network).gif").path, fps = 30)
end