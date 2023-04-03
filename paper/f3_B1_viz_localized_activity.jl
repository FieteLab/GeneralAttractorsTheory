
"""
Visualize a bump of activity on the neural lattice
"""


include("settings.jl")


duration = 500
still = 15
network = "torus"

for network in networks
    network == "sphere" && continue
    network != "plane" && continue
    print(hLine("Doing $network"; style ="bold blue"))

    φ = embeddings[network]

    # simulation params
    x₀ = if network == "sphere"
        [1, 0, 0]
    elseif network == "cylinder"
        [3, 0] 
    elseif network ∉ ("line", "ring")
        [-9.9, -9.9]
    elseif network == "ring"
        [3,]
    else
        [0,]
    end

    # plots params
    can = make_single_can(network; k = LocalGlobalKernel(α = 1.5, σ = 250.0))

    # run simulation
    h, _ = simulate_constant_traj_random_init(can, duration, dt, still, τ, b₀; x₀=x₀)

    activity = reshape(h.S[:, 1, end], can.n)

    # plot
    fig = if can.d == 2
        Δ = 2
        
        n_idx = 200

        coords = by_column(φ, can.X[:, 1:Δ:end])

        w_x = range(can.X[1,1], can.X[1,end]; length=can.n[1])
        w_y = range(can.X[2,1], can.X[2,end]; length=can.n[2])
        activity = reshape(h.S[:, 1, end], can.n)

        main_plot = contourf(
            w_x, w_y,
            activity, title = "Activity ($network)",
            aspect_ratio = :equal,
            color=:bwr,
            linewidth = 0.25,
            msc=:black,
            lc=:black,
            grid = false,
            levels=3,
            )

        top_plot = plot(w_x, sum(activity; dims=1)[1, :], label=nothing; 
                lw=3, color=:black, grid=false,
                xlabel="θ", ylabel="activity"
            )

            
        side_plot = plot(
                sum(activity; dims=2), w_y, label=nothing; 
                lw=3, color=:black, grid=false,
                xlabel="activity", ylabel="θ"
        )
        

        plot(
            plot(; plot_remove_axes...), 
            top_plot, 
            side_plot, 
            main_plot, 
            # grid(2, 2, heights=[.2, .8], widts=[.2, .8]),
            size=(1000, 1000)
        )
        

    else
        fig = plot(
            can.X[1, :],  
            activity,
            lw=2, color=:black, 
            fillrange = 0,
            fillalpha = 0.25,
            fillcolor = :black,
            label=nothing, xlabel="θ", ylabel="s(θ)",
            title = " $(network))",
        )

    end


    display(fig)
    save_plot(supervisor, fig, "f3_B_loc_actvitity_$(network)");
end