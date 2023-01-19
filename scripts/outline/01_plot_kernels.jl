include("settings.jl")
import Statistics: mean


parameters_range = kernels_parameters_range["torus"]

σs = Dict(
    :mexican_hat => 50, 
    :DoE => 100,
    :local_global => 50,
    :constant => 5,
)

plts =[]
for (name, ktype) in pairs(kernels)
    P = plot()
    # plot a few with random params
    for i in 1:25
        params = Dict(
            k => rand(v) for (k, v) in parameters_range[name]
        )
        plot!(
            ktype(; params...); σ=σs[name], lw=1, color=:gray, alpha=.5)
    end

    hline!([0], color=:blue, lw=2, ls=:dash, alpha=.2, label=nothing)
    vline!([0], color=:blue, lw=2, ls=:dash, alpha=.2, label=nothing)

    # plot mean values
    params = Dict(
        k => mean(v) for (k, v) in parameters_range[name]
    )
    plot!(
        ktype(; params...); σ=σs[name], lw=4, color=:black, alpha=1)
    push!(plts, P)
end

fig = plot(plts..., layout=(2,length(kernels)), size=(900, 500), 
        # ylim=[-3, 3], xlim=[-2.5, 2.5]
)

save_plot(supervisor, fig, "01_kernels")


fig