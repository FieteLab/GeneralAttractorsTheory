using Plots
using GeneralAttractors.Kernels
using GeneralAttractors.ProjectSupervisor

include("settings.jl")


"""
Plot localized activity peaks for the torus only.
"""

# can params
can_type = :single
kernel_name = :constant
can_name = "torus"

m = 48
φ = embeddings[can_name]
can_maker = network_makers[can_name]

# simulation params
duration = 25
still = 15

# make CAN
params = Dict(
    k => mean(v) for (k, v) in kernels_parameters_range[can_name][kernel_name]
)
can = can_maker(can_type; n=m, k=kernels[kernel_name](; params...))


plots = []
for i in 1:10
    x₀ = if i == 1
        [3.14, 3.14]
    else
        nothing
    end
    # x₀ = nothing

    # run simulation
    h, _ = simulate_constant_traj_random_init(can, duration, dt, still, τ, b₀; x₀=x₀)

    # plot

    push!(plots, 
        scatter3d(
            can.X[1, :], can.X[2, :], h.S[:, 1, end],
            marker_z = h.S[:, 1, end], msa=0, msw=0,
            colorbar = i == 5,
            label=nothing, 
            # clims=(-0.2, .8),
            showaxis = false,
            axis=nothing,
        )
    )
end

fig = plot(plots..., layout=(2, 5), size=(1200, 240))

save_plot(supervisor, fig, "03b_torus_localized_activity_$(kernel_name)")