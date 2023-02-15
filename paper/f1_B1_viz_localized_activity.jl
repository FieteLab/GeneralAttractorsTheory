
include("settings.jl")


duration = 25
still = 15

for can_name in networks
    print(hLine("Doing $can_name"; style ="bold blue"))
    
    φ = embeddings[can_name]
    can_maker = network_makers[can_name]

    # simulation params

    x₀ = if can_name == "sphere"
        [1, 0, 0]
    elseif can_name == "cylinder"
        [2π, -1] 
    elseif name ∉ ("line", "ring")
        [1, 0]
    else
        [1,]
    end

    # plots params
    can = can_maker(can_type;)

    # run simulation
    h, _ = simulate_constant_traj_random_init(can, duration, dt, still, τ, b₀; x₀=x₀)

    # plot
    fig = if can.d == 2
        Δ = 2
        activity = reshape(h.S[:, 1, end], can.n)
        coords = by_column(φ, can.X[:, 1:Δ:end])

        w_x = range(can.X[1,1], can.X[1,end]; length=can.n[1])
        w_y = range(can.X[2,1], can.X[2,end]; length=can.n[2])


        x = can.X[1,:]

        # push!(plots, heatmap(activity, title=name))
        # @info "cacca" can.X activity h.S[:, 1, end]
        # plt = scatter3d(
        #     can.X[1, :], can.X[2, :], h.S[:, 1, end],
        #     marker_z = h.S[:, 1, end], msa=0, msw=0
        # )
        error("make contourf")
    else

    end



    save_plot(supervisor, _fig, "f1_B_loc_actvitity_$(can_name)");
    break
end