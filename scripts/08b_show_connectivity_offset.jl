include("settings.jl")


plots = []

neuron_n = 1000

for δ in (0.0, 0.2)
    can = network_makers["torus"](:default; n=48, offset_size = δ, k = LocalGlobalKernel(α = 2.5, σ = 1.0))
    n = δ == 0.0 ? 1 : 4
    for i in 1:n
        
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
        push!(plots, plt)
        
    end
end

fig = plot(plots..., layout=(1, 5), size=(1000, 800))

save_plot(supervisor, fig, "08b_connectivity_offset")