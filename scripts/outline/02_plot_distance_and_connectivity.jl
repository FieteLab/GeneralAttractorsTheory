using GeneralAttractors
using Distances
import GeneralAttractors: plot_distance_function, SphericalDistance, MobiusEuclidean

d = PeriodicEuclidean([2π])
p = plot_distance_function(d, plot_title = "Ring distance")
save_plot(supervisor, p, "02_distance_ring")

d = PeriodicEuclidean([2π, 2π])  # distance function over a torus manifold
p = plot_distance_function(d, plot_title = "Torus distance")
save_plot(supervisor, p, "02_distance_torus")

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
