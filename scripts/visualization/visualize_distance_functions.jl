using GeneralAttractors
using Distances


# d = PeriodicEuclidean([2π])
# plot_distance_function(d, plot_title = "Ring distance")

# d = PeriodicEuclidean([2π, 2π])  # distance function over a torus manifold
# plot_distance_function(d, plot_title = "Torus distance")

# d = PeriodicEuclidean([2π, Inf])
# plot_distance_function(d, plot_title = "Cylinder distance", layout = (4, 1))

d = MobiusEuclidean()
plot_distance_function(d, plot_title = "Mobius distance", layout = (1, 4), colorbar=nothing)

# d = SphericalAngle()
# plot_distance_function(
#     d,
#     plot_title = "Spherical (angle) distance",
#     layout = (2, 2),
#     size = (800, 600),
# )

# import GeneralAttractors: SphericalDistance
# d = SphericalDistance()
# plot_distance_function(
#     d,
#     plot_title = "Spherical (angle) distance",
#     layout = (2, 2),
#     size = (800, 600),
# )
