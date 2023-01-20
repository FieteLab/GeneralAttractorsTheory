using Plots
using GeneralAttractors.Kernels
using GeneralAttractors.ProjectSupervisor

include("settings.jl")

# can params
can_type = :single
can_name = "torus"  # make sure it matches the can maker
m = 64
φ = embeddings[can_name]
can_maker = network_makers[can_name]

# simulation params
duration = 125
still = 100
x₀ = [1, 0]

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
        can = can_maker(can_type; m=m, k=ktype(; params...))

        # run simulation
        h, _ = simulate_constant_traj_random_init(can, duration, dt, still, τ, b₀; x₀=x₀)

        # plot
        Δ = 5
        activity = reshape(h.S[:, 1, end], can.n)
        coords = by_column(φ, can.X[:, 1:Δ:end])
        # push!(plots, heatmap(activity, title=name))
        push!(plots, 
        
            scatter3d(eachrow(coords)...;
                msa=0, msw=0, alpha=1,
                ms=10,
                xlim=(-1.25, 1.25), ylim=(-1.25, 1.25), zlim=(-1.25, 1.25),
                marker_z=h.S[1:Δ:end, 1, end],
                camera=(45, 50),
                label=nothing,
                title=[name],
                showaxis = false,
                axis=nothing,
            )
        )
    end
end

fig = plot(plots..., 
        clims=(-0.2, 0.8),
        layout=(length(kernels), n_samples_per_kernel), 
        size=(1400, 1400),
        aspect_ratio=:equal
)

save_plot(supervisor, fig, "03_actvitity");

# fig
nothing