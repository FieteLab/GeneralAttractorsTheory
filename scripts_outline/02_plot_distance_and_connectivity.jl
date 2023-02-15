using GeneralAttractors
using Distances
import GeneralAttractors: plot_distance_function, SphericalDistance, MobiusEuclidean


# ----------------------------- plot connectivity ---------------------------- #
plots = []

for network in networks
    can = network_makers[network](:single)

    W = reshape(can.W[1,:], can.n)' |> Matrix
    w_x = range(can.X[1,1], can.X[1,end]; length=can.n[1])
    w_y = range(can.X[2,1], can.X[2,end]; length=can.n[2])


    x = can.X[1,:]

    plt = contourf(
        w_x, w_y,
        W, title = "Connectivity matrix ($network)",
        aspect_ratio = :equal,
        color=:bwr,
        linewidth = 0.25,
        msc=:black,
        lc=:black,
        grid = false,
        )

    scatter!(plt, [x[1]], [x[2]], color=:green, ms=10, msw=0, msa=0, label=nothing)
    push!(plots, plt)
end

fig = plot(plots..., size=(1000, 1000))
save_plot(supervisor, fig, "02_connectivity")
display(fig)


# ------------------------------- plot distance ------------------------------ #

d = PeriodicEuclidean([2π])
p = plot_distance_function(d, plot_title = "Ring distance")
save_plot(supervisor, p, "02_distance_ring")

d = PeriodicEuclidean([2π, 2π])  # distance function over a torus manifold
p = plot_distance_function(d, plot_title = "Torus distance")
save_plot(supervisor, p, "02_distance_torus")


d = PeriodicEuclidean([Inf, Inf])  
p = plot_distance_function(d, plot_title = "Plane distance")
save_plot(supervisor, p, "02_distance_plane")

d = PeriodicEuclidean([2π, Inf])
p = plot_distance_function(d, plot_title = "Cylinder distance", layout = (4, 1))
save_plot(supervisor, p, "02_distance_cylinder")

d = MobiusEuclidean()
p = plot_distance_function(
    d,
    plot_title = "Mobius distance",
    layout = (1, 4),
    colorbar = nothing,
)
save_plot(supervisor, p, "02_distance_mobius")

# d = SphericalDistance()
# plot_distance_function(
#     d,
#     plot_title = "Spherical (angle) distance",
#     layout = (2, 2),
#     size = (800, 600),
# )
