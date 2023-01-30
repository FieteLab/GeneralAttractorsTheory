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
temporal_downsampling = 50
fps = 20


for network in ("sphere", "mobius")

    savepath = supervisor.projectdir / "plots" / "path_int_$(network).gif"
    exists(savepath) && continue

    # make network to get neurons lattice coordinates
    can = network_makers[network](:single; n=48)
    w_x = range(can.X[1,1], can.X[1,end]; length=can.n[1])
    w_y = range(can.X[2,1], can.X[2,end]; length=can.n[2])

    d = can.name == "sphere" ? 3 : 2
    coord3d = d == 2 ? by_column(embeddings[network], can.X) : can.X

    # --------------------------------- get data --------------------------------- #
    # load 3d iso embedding model and embedded data
    set_datadir(supervisor,  (supervisor.datadir / "mfld_top").path)
    embs, M = ProjectSupervisor.fetch(
        supervisor; 
        tag="d3_embeddings", 
        can=network, 
        dim=3
    )[2]
    pca = embs[:pca]
    iso = embs[:iso]
    @assert size(M, 1) == 3 size(M)


    # load path integration data
    set_datadir(supervisor,  ((supervisor.datadir-1) / "path_int").path)
    filters = Dict{Symbol, Any}(
        :tag => tag,
        :network => network,
        :duration => 2500,
    )

    data = ProjectSupervisor.fetch(supervisor; filters...)[2][1]
    history = data["h"]
    trajectory = data["trajectory"]
    S = sum(history.S; dims=2)[:, 1, :]
    S_embedd = predict(iso, predict(pca, S))
    X = history.x_M_decoded


    @info "Got data for $network" S M iso

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

            frame_content = path_int_gif_frame(trajectory, X, S, can, w_x, w_y, M, S_embedd, fnum, skipframes)

            frame(anim)
            update!(job)
        end
    end

    gif(anim, (savepath).path, fps = fps)
end