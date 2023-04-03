include("settings.jl")

plots = []
network = "sphere"
neuron_n = 1600

# for δ in (0.0, 0.2)
δ = .25
can = network_makers[network](:default; n=48, offset_size = δ)
n = δ == 0.0 ? 1 : 6
for i in 1:n
    if network != "sphere"
        W = reshape(can.Ws[i][neuron_n,:], can.n)' |> Matrix
        w_x = range(can.X[1,1], can.X[1,end]; length=can.n[1])
        w_y = range(can.X[2,1], can.X[2,end]; length=can.n[2])


        x = can.X[:, neuron_n]

        plt = contourf(
            w_x, w_y,
            W, title = "Connectivity matrix δ=$δ",
            aspect_ratio = :equal,
            color=:bwr,
            linewidth = 0.25,
            msc=:black,
            lc=:black,
            grid = false,
            colorbar = false,
            # clims=(-1, 0)
            )

        scatter!(plt, [x[1]], [x[2]], color=:green, ms=2, msw=0, msa=0, label=nothing)
    else
        plt = scatter(
            eachrow(can.X)...,
            marker_z=can.Ws[i][neuron_n,:],
            aspect_ratio = :equal,
            color=:bwr,
            linewidth = 0.25,
            msa=0, msw=0,
            ms=4,
            lc=:black,
            grid = false,
            colorbar = false,
            legend=false,
            clims=(-.05, 0),
            xlim=[-1.1, 1.1],
            ylim=[-1.1, 1.1],
            zlim=[-1.1, 1.1],
        )

        scatter!(
            plt,
            [can.X[1, neuron_n]],
            [can.X[2, neuron_n]],
            [can.X[3, neuron_n]],
            color=:green,
            ms=10,
            msw=0,
            msa=0,
            label=nothing,
        )
    end
    push!(plots, plt)
end
# end

fig = plot(plots..., layout=(3, 2), size=(1000,1000),)
display(fig)
save_plot(supervisor, fig, "f5_C1_offset_connectivity_$(network)")