
function animate_3d_scatter(
        X::Matrix, 
            savefld::String="test",
            savename::String="test";
            downsample=1, alpha=.5, ms=5,
            roll_fn::Function = (t)->t*2,
            azimuth_fn::Function = (t)->50,
            nframes=100, fps=15, kwargs...
            )
    @assert size(X, 1) == 3 "X size wrong: $(size(X)) should have second dimension == 3"

    # h = maximum(abs.(X))
    p = scatter(
        X[1, 1:downsample:end], X[2, 1:downsample:end], X[3, 1:downsample:end], 
        alpha=alpha, label=nothing, ms=ms,
        color=:black,  msw=0.0
        # ylim=[-h, h], xlim=[-h, h], zlim=[-h, h]
        ; kwargs...
    )
    anim = Animation()

    pbar = ProgressBar()
    job = addjob!(pbar, description="Making animation",  N=nframes)
    Progress.with(pbar) do
        for i in 1:nframes
            plot!(camera=(roll_fn(i), azimuth_fn(i)))
            frame(anim)
            update!(job)
        end
    end

    gif(anim, savepath(savefld, savename, "gif"), fps=fps)
end