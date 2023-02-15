include("settings.jl")
import Statistics: mean


parameters_range = kernels_parameters_range["torus"]


this_plot_params_ranges = Dict(
    :mexican_hat => Dict(
        :α => 1:δ:7,
        :σ => 5:δ:20,
    ),
    :DoE => Dict(
        :a => 1:δ:3,
        :λ => 5:δ:10,
    ),
    :local_global => Dict(
        :α => 0.5:δ:1.5,
        :σ => 5:δ:10,
    ),
    :constant => Dict(
        :σ => 1.0:δ:10.0,
        :β_minus => -1.8:δ:-.1
)
)

plts =[]
for (name, ktype) in pairs(kernels)
    P = plot()
    # plot a few with random params
    for i in 1:25
        params = Dict(
            k => rand(v) for (k, v) in this_plot_params_ranges[name]
        )
        plot!(
            ktype(; params...); σ=15, lw=1, color=:gray, alpha=.5)
    end

    hline!([0], color=:blue, lw=2, ls=:dash, alpha=.2, label=nothing)
    vline!([0], color=:blue, lw=2, ls=:dash, alpha=.2, label=nothing)

    # plot mean values
    params = Dict(
        k => mean(v) for (k, v) in this_plot_params_ranges[name]
    )
    plot!(
        ktype(; params...); σ=15, lw=4, color=:black, alpha=1,
        ylim=[-1.75, 0.25], yticks=[-1, 0], 
        xticks = [-10, 0, 10],
        title=name, label=nothing
        )
    push!(plts, P)
end

fig = plot(plts..., layout=(2,length(kernels)), size=(900, 500),)

save_plot(supervisor, fig, "01_kernels")


fig