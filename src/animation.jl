using Plots
using Term.Progress

function animate_simulation_data(can, traj, hist, X, φ, savepath; dt=0.5, frames_Δ=25, neurons_Δ = 11)
    coord3d = by_column(φ, can.X)
    nframes = size(traj.X, 1)
    skipframes = traj.still / dt |> floor |> Int

    anim = Animation()

    pbar = ProgressBar()
    Progress.with(pbar) do
        job = addjob!(pbar, description = "Making animation", N = (Int ∘ ceil)((nframes - skipframes) / frames_Δ))
        for fnum in 1:frames_Δ:nframes
            fnum * dt < traj.still && continue
            # fnum % hist.average_over == 0 && (hist_frame += 1)
            hist_frame = (fnum-skipframes)  ÷ hist.average_over + 1
            # plot trajectory
            try
                p1 = plot(traj, fnum)
                plot!(
                    eachcol(X[1:fnum, :])...,
                    lw = 6, color=:red, alpha=.5,
                    label = "decoded",
                    )

                # plot neurons activation
                s = sum(hist.S[:, :, hist_frame], dims=2)[1:neurons_Δ:end]
                
                # p2 = scatter3d(
                #     coord3d[1, 1:neurons_Δ:end], coord3d[2, 1:neurons_Δ:end], coord3d[3, 1:neurons_Δ:end],
                #     marker_z = s,
                #     ms = 5,  msw = 0.0,
                #     alpha=0.5,
                #     xlim = [minimum(coord3d[1,:])-.5, maximum(coord3d[1,:])+.5], 
                #     ylim = [minimum(coord3d[2,:])-.5, maximum(coord3d[2,:])+.5], 
                #     zlim = [minimum(coord3d[3,:])-.5, maximum(coord3d[3,:])+.5], 
                #     size = (800, 800),
                #     colorbar=false,
                #     )

                p2 = scatter(
                    eachrow(can.X[:, 1:neurons_Δ:end])...,
                    marker_z=s,
                    )

    
                # plot(simulation, time[framen], framen, x, v, X̄, φ)
                plot(p1, p2, layout = (1, 2), size = (1200, 800))
                frame(anim)
            catch e
                break
            end
            update!(job)
        end
        
    end

    gif(anim, savepath, fps = 30)
end