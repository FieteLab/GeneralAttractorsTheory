using Plots
using GeneralAttractors.Kernels
using GeneralAttractors.ProjectSupervisor

include("settings.jl")

# can params
can_type = :single

for can_name in networks
    print(hLine("Doing $can_name"; style ="bold blue"))
    
    m = 48
    φ = embeddings[can_name]
    can_maker = network_makers[can_name]

    # simulation params
    duration = 25
    still = 15
    x₀ = if can_name == "sphere"
        [1, 0, 0]
    elseif can_name == "cylinder"
        [2π, -1] 
    else
        [1, 0]
    end
    # x₀ = nothing

    # plots params
    n_samples_per_kernel = 2
    plots = []

    for (name, ktype) in pairs(kernels)
        @info "doing $name"
        ranges = kernels_parameters_range[can_name][name]

        for i in 1:n_samples_per_kernel
            # make CAN
            params = Dict(
                k => i == 1 ? minimum(v) : maximum(v) for (k, v) in ranges
            )
            can = can_maker(can_type; n=m, k=ktype(; params...))

            # run simulation
            h, _ = simulate_constant_traj_random_init(can, duration, dt, still, τ, b₀; x₀=x₀)

            # plot
            Δ = 2
            activity = reshape(h.S[:, 1, end], can.n)
            coords = by_column(φ, can.X[:, 1:Δ:end])
  
            w_x = range(can.X[1,1], can.X[1,end]; length=can.n[1])
            w_y = range(can.X[2,1], can.X[2,end]; length=can.n[2])
        
        
            x = can.X[1,:]
        
            # push!(plots, heatmap(activity, title=name))
            # @info "cacca" can.X activity h.S[:, 1, end]
            plt = scatter3d(
                can.X[1, :], can.X[2, :], h.S[:, 1, end],
                marker_z = h.S[:, 1, end], msa=0, msw=0
            )

            # push!(plots, 
            #     scatter3d(eachrow(coords)...;
            #         msa=0, msw=0, alpha=1,
            #         ms=8,
            #         xlim=(-1.25, 1.25), ylim=(-1.25, 1.25), zlim=(-1.25, 1.25),
            #         marker_z=h.S[1:Δ:end, 1, end],
            #         camera=(45, 50),
            #         label=nothing,
            #         title=[name],
            #         showaxis = false,
            #         axis=nothing,
            #     )
            # )
        end
    end

    _fig = plot(plots..., 
            clims=(-0.2, 0.8),
            layout=(n_samples_per_kernel, length(kernels)), 
            size=(1400, 1000),
            aspect_ratio=:equal
    )
    # save_plot(supervisor, _fig, "03_actvitity_$(can_name)");
    display(_fig)

    break
end

# nothing