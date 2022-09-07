using GeneralAttractors
using Distances


d = PeriodicEuclidean([2π])
plot_distance_function(d, plot_title="Ring distance")

d = PeriodicEuclidean([2π, 2π])  # distance function over a torus manifold
plot_distance_function(d, plot_title="Torus distance")

d = PeriodicEuclidean([2π, Inf])  
plot_distance_function(d, plot_title="Cylinder distance", layout=(2, 2))