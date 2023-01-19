using Plots
using GeneralAttractors.Kernels
using GeneralAttractors.ProjectSupervisor

sup = Supervisor("GeneralAttractorsTheory")

kernels = (;
   :mexican_hat => MexicanHatKernel, 
   :DoE => DiffOfExpKernel,
   :local_global => LocalGlobalKernel,
   :constant => ConstantKernel,
)
parameters_range = Dict(
    :mexican_hat => Dict(
        :α => 0.1:0.01:0.5,
        :σ => 0.1:0.01:0.5,
    ),
    :DoE => Dict(
        :a => 1.25:0.01:3,
        :λ => 0.5:0.01:2.75,
    ),
    :local_global => Dict(
        :α => 0.5:.01:1.5,
        :σ => 0.1:.01:1.5,
    ),
    :constant => Dict(
        :σ => 0.5:.01:1.5,
    ),
)

plts =[]
for (name, ktype) in pairs(kernels)
    P = plot()
    # plot a few with random params
    for i in 1:25
        params = Dict(
            k => rand(v) for (k, v) in parameters_range[name]
        )
        plot!(ktype(; params...); lw=1, color=:gray, alpha=.5)
    end

    hline!([0], color=:black, lw=2, ls=:dash, alpha=.2, label=nothing)

    # plot default values
    plot!(ktype(); title=name, lw=4, color=:black)
    push!(plts, P)
end

fig = plot(plts..., layout=(2,length(kernels)), size=(900, 300), 
        ylim=[-3, 3], xlim=[-2.5, 2.5]
)

save_plot(sup, fig, "01_kernels")